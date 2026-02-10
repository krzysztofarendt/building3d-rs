use anyhow::Result;

use crate::sim::energy::config::ThermalConfig;
use crate::sim::energy::hvac::HvacIdealLoads;
use crate::sim::materials::MaterialLibrary;
use crate::{Building, UID};

use super::ThermalNetwork;
use super::multizone::MultiZoneStepResult;
use super::solve::solve_with_fixed_nodes;
use crate::sim::energy::boundary::ThermalBoundaries;
use crate::sim::index::SurfaceIndex;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Control {
    None,
    Heat,
    Cool,
}

/// Multi-zone air-node model with a single 2R1C envelope node per zone.
///
/// Each zone has:
/// - an air temperature node `T_air`,
/// - an envelope temperature node `T_env` representing aggregated wall mass.
///
/// The steady-state exterior conductance `K_env = Σ(U*A)` is split equally:
/// - `K_air_env = 2*K_env` between air and envelope,
/// - `K_env_out = 2*K_env` between envelope and outdoor.
///
/// This preserves the same steady-state heat loss to outside while introducing
/// thermal mass dynamics through `C_env`.
#[derive(Debug, Clone)]
pub struct MultiZoneEnvelopeRcModel {
    zone_uids: Vec<UID>,
    zone_names: Vec<String>,
    air_temperatures_c: Vec<f64>,
    env_temperatures_c: Vec<f64>,
    air_capacity_j_per_k: Vec<f64>,
    env_capacity_j_per_k: Vec<f64>,
    exterior_k_w_per_k: Vec<f64>,
    infiltration_k_w_per_k: Vec<f64>,
    neighbors_air: Vec<Vec<(usize, f64)>>,
    sum_interzone_k: Vec<f64>,
}

impl MultiZoneEnvelopeRcModel {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        building: &Building,
        network: &ThermalNetwork,
        surface_index: &SurfaceIndex,
        boundaries: &ThermalBoundaries,
        infiltration_ach: f64,
        air_capacity_j_per_m3_k: f64,
        default_envelope_capacity_j_per_m2_k: f64,
        material_library: Option<&MaterialLibrary>,
        initial_temp_c: f64,
    ) -> Self {
        let cap_area_for_surface = |surface: &crate::sim::index::SurfaceRef| -> f64 {
            material_library
                .and_then(|lib| lib.lookup(&surface.path))
                .and_then(|m| m.thermal.as_ref())
                .map(|t| t.thermal_capacity)
                .unwrap_or(default_envelope_capacity_j_per_m2_k)
                .max(0.0)
        };

        Self::new_impl(
            building,
            network,
            surface_index,
            boundaries,
            infiltration_ach,
            air_capacity_j_per_m3_k,
            cap_area_for_surface,
            initial_temp_c,
        )
    }

    /// Like [`Self::new`], but resolves envelope capacity per surface using the shared
    /// [`ThermalConfig`] precedence rules (UID override > exact path construction >
    /// MaterialLibrary thermal capacity > pattern construction > default).
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_thermal_config(
        building: &Building,
        network: &ThermalNetwork,
        surface_index: &SurfaceIndex,
        boundaries: &ThermalBoundaries,
        thermal_config: &ThermalConfig,
        infiltration_ach: f64,
        air_capacity_j_per_m3_k: f64,
        default_envelope_capacity_j_per_m2_k: f64,
        initial_temp_c: f64,
    ) -> Self {
        let cap_area_for_surface = |surface: &crate::sim::index::SurfaceRef| -> f64 {
            thermal_config.resolve_envelope_capacity_j_per_m2_k(
                Some(&surface.polygon_uid),
                &surface.path,
                default_envelope_capacity_j_per_m2_k,
            )
        };

        Self::new_impl(
            building,
            network,
            surface_index,
            boundaries,
            infiltration_ach,
            air_capacity_j_per_m3_k,
            cap_area_for_surface,
            initial_temp_c,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn new_impl(
        building: &Building,
        network: &ThermalNetwork,
        surface_index: &SurfaceIndex,
        boundaries: &ThermalBoundaries,
        infiltration_ach: f64,
        air_capacity_j_per_m3_k: f64,
        cap_area_for_surface: impl Fn(&crate::sim::index::SurfaceRef) -> f64,
        initial_temp_c: f64,
    ) -> Self {
        // Stable zone ordering: building.zones() is name-sorted by convention.
        let zones = building.zones();
        let n = zones.len();

        let mut zone_uids = Vec::with_capacity(n);
        let mut zone_names = Vec::with_capacity(n);
        let mut zone_volumes_m3 = Vec::with_capacity(n);
        for z in zones {
            zone_uids.push(z.uid.clone());
            zone_names.push(z.name.clone());
            zone_volumes_m3.push(z.volume());
        }

        let air_temperatures_c = vec![initial_temp_c; n];
        let env_temperatures_c = vec![initial_temp_c; n];

        let air_capacity_j_per_k: Vec<f64> = zone_volumes_m3
            .iter()
            .map(|v| v * air_capacity_j_per_m3_k.max(0.0))
            .collect();

        // Exterior conductance per zone (Σ U*A).
        let exterior_k_w_per_k: Vec<f64> = zone_uids
            .iter()
            .map(|uid| network.exterior_conductance_by_zone_w_per_k(uid))
            .collect();

        // Infiltration conductance: rho * cp * V * ACH / 3600.
        let infiltration_k_w_per_k: Vec<f64> = zone_volumes_m3
            .iter()
            .map(|v| 1.2 * 1005.0 * v * infiltration_ach / 3600.0)
            .collect();

        // Envelope capacity per zone: Σ (C_area(path) * area).
        let mut env_capacity_j_per_k = vec![0.0; n];
        let uid_to_idx: std::collections::HashMap<&str, usize> = zone_uids
            .iter()
            .enumerate()
            .map(|(i, uid)| (uid.as_str(), i))
            .collect();
        for surface in &surface_index.surfaces {
            if !boundaries.is_exterior(&surface.polygon_uid) {
                continue;
            }
            let Some(&i) = uid_to_idx.get(surface.zone_uid.as_str()) else {
                continue;
            };

            let cap_area = cap_area_for_surface(surface);

            env_capacity_j_per_k[i] += cap_area * surface.area_m2;
        }

        // Inter-zone adjacency list for air nodes.
        let mut neighbors_air: Vec<Vec<(usize, f64)>> = vec![vec![]; n];
        let mut sum_interzone_k = vec![0.0; n];
        let uid_to_idx2: std::collections::HashMap<&str, usize> = zone_uids
            .iter()
            .enumerate()
            .map(|(i, uid)| (uid.as_str(), i))
            .collect();

        for e in network.interzone_conductances() {
            let Some(&i) = uid_to_idx2.get(e.zone_a.as_str()) else {
                continue;
            };
            let Some(&j) = uid_to_idx2.get(e.zone_b.as_str()) else {
                continue;
            };
            if i == j {
                continue;
            }
            neighbors_air[i].push((j, e.conductance_w_per_k));
            neighbors_air[j].push((i, e.conductance_w_per_k));
            sum_interzone_k[i] += e.conductance_w_per_k;
            sum_interzone_k[j] += e.conductance_w_per_k;
        }

        Self {
            zone_uids,
            zone_names,
            air_temperatures_c,
            env_temperatures_c,
            air_capacity_j_per_k,
            env_capacity_j_per_k,
            exterior_k_w_per_k,
            infiltration_k_w_per_k,
            neighbors_air,
            sum_interzone_k,
        }
    }

    pub fn zone_uids(&self) -> &[UID] {
        &self.zone_uids
    }

    pub fn zone_names(&self) -> &[String] {
        &self.zone_names
    }

    pub fn air_temperatures_c(&self) -> &[f64] {
        &self.air_temperatures_c
    }

    pub fn envelope_temperatures_c(&self) -> &[f64] {
        &self.env_temperatures_c
    }

    /// Advances the model by one timestep.
    ///
    /// - `outdoor_temp_c`: boundary temperature for exterior + infiltration conductances
    /// - `air_gains_w`: per-zone gains applied to the *air* node (W), same ordering as `zone_uids()`
    /// - `envelope_gains_w`: per-zone gains applied to the *envelope* node (W), same ordering as `zone_uids()`
    /// - `hvac`: ideal-loads controller (setpoints)
    /// - `dt_s`: timestep in seconds
    pub fn step(
        &mut self,
        outdoor_temp_c: f64,
        air_gains_w: &[f64],
        envelope_gains_w: &[f64],
        hvac: &HvacIdealLoads,
        dt_s: f64,
    ) -> Result<MultiZoneStepResult> {
        let n = self.zone_uids.len();
        anyhow::ensure!(air_gains_w.len() == n, "air gains length mismatch");
        anyhow::ensure!(
            envelope_gains_w.len() == n,
            "envelope gains length mismatch"
        );
        anyhow::ensure!(dt_s > 0.0, "dt must be positive");

        // Unknown ordering: [air_0..air_{n-1}, env_0..env_{n-1}]
        let nn = 2 * n;
        let mut diag = vec![0.0; nn];
        let mut rhs_base = vec![0.0; nn];
        let mut neighbors: Vec<Vec<(usize, f64)>> = vec![vec![]; nn];

        // Build air rows.
        for (i, &air_gain_w) in air_gains_w.iter().enumerate() {
            let c_air = self.air_capacity_j_per_k[i].max(0.0);
            let c_over_dt = if c_air > 0.0 { c_air / dt_s } else { 0.0 };
            let k_inf = self.infiltration_k_w_per_k[i];
            let k_env = self.exterior_k_w_per_k[i].max(0.0);
            let k_air_env = 2.0 * k_env;

            let idx_air = i;
            let idx_env = n + i;

            diag[idx_air] = c_over_dt + k_inf + k_air_env + self.sum_interzone_k[i];
            rhs_base[idx_air] =
                c_over_dt * self.air_temperatures_c[i] + k_inf * outdoor_temp_c + air_gain_w;

            // Air ↔ envelope coupling.
            if k_air_env > 0.0 {
                neighbors[idx_air].push((idx_env, k_air_env));
            }

            // Air ↔ air coupling (inter-zone).
            for &(j, k) in &self.neighbors_air[i] {
                if k > 0.0 {
                    neighbors[idx_air].push((j, k));
                }
            }
        }

        // Build envelope rows.
        for (i, &envelope_gain_w) in envelope_gains_w.iter().enumerate() {
            let c_env = self.env_capacity_j_per_k[i].max(0.0);
            let c_over_dt = if c_env > 0.0 { c_env / dt_s } else { 0.0 };
            let k_env = self.exterior_k_w_per_k[i].max(0.0);
            let k_air_env = 2.0 * k_env;
            let k_env_out = 2.0 * k_env;

            let idx_env = n + i;
            let idx_air = i;

            diag[idx_env] = c_over_dt + k_air_env + k_env_out;
            rhs_base[idx_env] = c_over_dt * self.env_temperatures_c[i]
                + k_env_out * outdoor_temp_c
                + envelope_gain_w;

            if k_air_env > 0.0 {
                neighbors[idx_env].push((idx_air, k_air_env));
            }
        }

        // Guard against ill-posed air nodes (e.g., steady-state + zero conductances).
        // Envelope nodes may be degenerate (no exterior surfaces); those are fixed in the solver.
        for (i, &d) in diag[..n].iter().enumerate() {
            anyhow::ensure!(d > 0.0, "Ill-posed RC system: air diag[{i}] == 0");
        }

        let mut control = vec![Control::None; n];
        let mut setpoints = vec![0.0; n];

        // Iteratively determine which air nodes need setpoint control, accounting for coupling.
        for _ in 0..20 {
            let temps = solve_rc_with_controls(
                &diag,
                &rhs_base,
                &neighbors,
                &control,
                &setpoints,
                n,
                outdoor_temp_c,
            )?;

            let mut new_control = vec![Control::None; n];
            let mut new_setpoints = vec![0.0; n];

            for i in 0..n {
                let c_air = self.air_capacity_j_per_k[i].max(0.0);
                let c_over_dt = if c_air > 0.0 { c_air / dt_s } else { 0.0 };
                let k_inf = self.infiltration_k_w_per_k[i];
                let k_env = self.exterior_k_w_per_k[i].max(0.0);
                let k_air_env = 2.0 * k_env;
                let denom = c_over_dt + k_inf + k_air_env + self.sum_interzone_k[i];

                let idx_env = n + i;
                let mut neighbor_term = 0.0;
                for &(j, k) in &self.neighbors_air[i] {
                    neighbor_term += k * temps[j];
                }
                neighbor_term += k_air_env * temps[idx_env];

                let t_free = if denom > 1e-14 {
                    (c_over_dt * self.air_temperatures_c[i]
                        + k_inf * outdoor_temp_c
                        + air_gains_w[i]
                        + neighbor_term)
                        / denom
                } else {
                    temps[i]
                };

                if t_free < hvac.heating_setpoint {
                    new_control[i] = Control::Heat;
                    new_setpoints[i] = hvac.heating_setpoint;
                } else if t_free > hvac.cooling_setpoint {
                    new_control[i] = Control::Cool;
                    new_setpoints[i] = hvac.cooling_setpoint;
                }
            }

            if new_control == control && new_setpoints == setpoints {
                let final_temps = solve_rc_with_controls(
                    &diag,
                    &rhs_base,
                    &neighbors,
                    &control,
                    &setpoints,
                    n,
                    outdoor_temp_c,
                )?;

                let mut zone_heating_w = vec![0.0; n];
                let mut zone_cooling_w = vec![0.0; n];

                for i in 0..n {
                    if control[i] == Control::None {
                        continue;
                    }

                    let idx_air = i;
                    let idx_env = n + i;

                    let c_air = self.air_capacity_j_per_k[i].max(0.0);
                    let c_over_dt = if c_air > 0.0 { c_air / dt_s } else { 0.0 };
                    let k_inf = self.infiltration_k_w_per_k[i];
                    let k_env = self.exterior_k_w_per_k[i].max(0.0);
                    let k_air_env = 2.0 * k_env;

                    let mut coupling = 0.0;
                    for &(j, k) in &self.neighbors_air[i] {
                        coupling += k * (final_temps[idx_air] - final_temps[j]);
                    }
                    coupling += k_air_env * (final_temps[idx_air] - final_temps[idx_env]);

                    // Q_hvac balances the air node equation.
                    let q_hvac = c_over_dt * (final_temps[idx_air] - self.air_temperatures_c[i])
                        + k_inf * (final_temps[idx_air] - outdoor_temp_c)
                        + coupling
                        - air_gains_w[i];

                    if q_hvac > 0.0 {
                        zone_heating_w[i] = q_hvac;
                    } else {
                        zone_cooling_w[i] = -q_hvac;
                    }
                }

                self.air_temperatures_c[..n].copy_from_slice(&final_temps[..n]);
                self.env_temperatures_c[..n].copy_from_slice(&final_temps[n..(2 * n)]);

                return Ok(MultiZoneStepResult {
                    zone_uids: self.zone_uids.clone(),
                    zone_names: self.zone_names.clone(),
                    zone_temperatures_c: self.air_temperatures_c.clone(),
                    zone_heating_w,
                    zone_cooling_w,
                });
            }

            control = new_control;
            setpoints = new_setpoints;
        }

        anyhow::bail!("RC control iteration did not converge");
    }
}

fn solve_rc_with_controls(
    diag: &[f64],
    rhs_base: &[f64],
    neighbors: &[Vec<(usize, f64)>],
    control: &[Control],
    setpoints: &[f64],
    n_zones: usize,
    outdoor_temp_c: f64,
) -> Result<Vec<f64>> {
    let nn = diag.len();
    anyhow::ensure!(nn == 2 * n_zones, "dimension mismatch");

    let mut is_fixed = vec![false; nn];
    let mut fixed_temp = vec![0.0; nn];

    // Fixed air nodes (setpoints).
    for i in 0..n_zones {
        if control[i] != Control::None {
            is_fixed[i] = true;
            fixed_temp[i] = setpoints[i];
        }
    }

    // Degenerate envelope nodes (no capacity and no conductance) can exist if a zone has
    // no exterior surfaces. Fix them to the outdoor boundary to keep the system solvable.
    for i in 0..n_zones {
        let idx_env = n_zones + i;
        if diag[idx_env] <= 0.0 {
            is_fixed[idx_env] = true;
            fixed_temp[idx_env] = outdoor_temp_c;
        }
    }

    solve_with_fixed_nodes(diag, rhs_base, neighbors, &is_fixed, &fixed_temp)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sim::energy::boundary::ThermalBoundaries;
    use crate::sim::energy::config::ThermalConfig;
    use crate::sim::energy::network::ThermalNetwork;
    use crate::sim::energy::weather::WeatherData;
    use crate::sim::index::SurfaceIndex;
    use crate::{Building, Solid, Zone};

    #[test]
    fn test_envelope_rc_introduces_lag_vs_air_only() {
        // Single zone, pure envelope conductance, step outdoor from 20C to 0C.
        let s = Solid::from_box(3.0, 3.0, 3.0, None, "s").unwrap();
        let z = Zone::new("z", vec![s]).unwrap();
        let building = Building::new("b", vec![z]).unwrap();

        let index = SurfaceIndex::new(&building);
        let boundaries = ThermalBoundaries::classify(&building, &index);
        let mut cfg = ThermalConfig::new();
        cfg.infiltration_ach = 0.0;
        cfg.default_u_value = 1.0;
        cfg.indoor_temperature = 20.0;
        cfg.outdoor_temperature = 20.0;

        let network = ThermalNetwork::build(&building, &cfg, &index, &boundaries);
        let hvac = HvacIdealLoads::with_setpoints(-1e9, 1e9); // HVAC off

        let mut rc = MultiZoneEnvelopeRcModel::new(
            &building,
            &network,
            &index,
            &boundaries,
            0.0,
            0.0,
            200_000.0,
            None,
            20.0,
        );

        // Outdoor temp drops; gains zero.
        let air_gains = vec![0.0];
        let env_gains = vec![0.0];
        let r0 = rc.step(0.0, &air_gains, &env_gains, &hvac, 3600.0).unwrap();
        let t_air_after_1 = r0.zone_temperatures_c[0];

        // With large envelope capacity, air shouldn't instantly hit the steady-state value.
        assert!(t_air_after_1 > 0.0);

        // Run several steps; should approach outdoor.
        for _ in 0..48 {
            rc.step(0.0, &air_gains, &env_gains, &hvac, 3600.0).unwrap();
        }
        assert!(rc.air_temperatures_c()[0] < 5.0);
    }

    #[test]
    fn test_rc_model_smoke_with_weather_synthetic() {
        let s = Solid::from_box(3.0, 3.0, 3.0, None, "s").unwrap();
        let z = Zone::new("z", vec![s]).unwrap();
        let building = Building::new("b", vec![z]).unwrap();
        let index = SurfaceIndex::new(&building);
        let boundaries = ThermalBoundaries::classify(&building, &index);
        let cfg = ThermalConfig::new();
        let network = ThermalNetwork::build(&building, &cfg, &index, &boundaries);

        let weather = WeatherData::synthetic("X", 52.0, 13.0, 10.0, 0.0);
        let mut model = MultiZoneEnvelopeRcModel::new(
            &building,
            &network,
            &index,
            &boundaries,
            cfg.infiltration_ach,
            cfg.thermal_capacity_j_per_m3_k,
            100_000.0,
            None,
            cfg.indoor_temperature,
        );
        let hvac = HvacIdealLoads::new();
        let air_gains = vec![0.0];
        let env_gains = vec![0.0];
        model
            .step(
                weather.records[0].dry_bulb_temperature,
                &air_gains,
                &env_gains,
                &hvac,
                3600.0,
            )
            .unwrap();
    }
}
