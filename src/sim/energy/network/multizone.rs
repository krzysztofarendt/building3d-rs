use anyhow::Result;

use super::ThermalNetwork;
use super::solve::solve_dense;
use crate::sim::energy::hvac::HvacIdealLoads;
use crate::{Building, UID};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Control {
    None,
    Heat,
    Cool,
}

/// Result of a single multi-zone timestep.
#[derive(Debug, Clone)]
pub struct MultiZoneStepResult {
    /// Zone UIDs in the same order as the arrays below.
    pub zone_uids: Vec<UID>,
    /// Zone names in the same order as the arrays below.
    pub zone_names: Vec<String>,
    /// Zone air temperatures after the step (°C).
    pub zone_temperatures_c: Vec<f64>,
    /// Thermal HVAC heating power per zone (W).
    pub zone_heating_w: Vec<f64>,
    /// Thermal HVAC cooling power per zone (W).
    pub zone_cooling_w: Vec<f64>,
}

/// Multi-zone air-node thermal model coupled by inter-zone conductances.
///
/// This is a “zone-air-only” model: each `Zone` is represented by a single air
/// temperature state, with coupling through partition conductances. This forms
/// the basis for future extensions (surface temperature nodes, radiant exchange,
/// airflow networks) without changing the coupling contracts.
#[derive(Debug, Clone)]
pub struct MultiZoneAirModel {
    zone_uids: Vec<UID>,
    zone_names: Vec<String>,
    temperatures_c: Vec<f64>,
    thermal_capacity_j_per_k: Vec<f64>,
    exterior_conductance_w_per_k: Vec<f64>,
    infiltration_conductance_w_per_k: Vec<f64>,
    neighbors: Vec<Vec<(usize, f64)>>,
    sum_interzone_k: Vec<f64>,
}

impl MultiZoneAirModel {
    /// Creates a model with initial zone temperatures equal to `initial_temp_c`.
    ///
    /// Thermal capacity is estimated from zone volume using a tuning factor of
    /// 50 kJ/(m³·K), consistent with the current energy module assumptions.
    pub fn new(
        building: &Building,
        network: &ThermalNetwork,
        infiltration_ach: f64,
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

        let temperatures_c = vec![initial_temp_c; n];

        // Tuning factor: 50 kJ/(m³·K) for medium-weight construction.
        let thermal_capacity_j_per_k: Vec<f64> =
            zone_volumes_m3.iter().map(|v| v * 50_000.0).collect();

        let exterior_conductance_w_per_k: Vec<f64> = zone_uids
            .iter()
            .map(|uid| network.exterior_conductance_by_zone_w_per_k(uid))
            .collect();

        // Infiltration conductance: rho * cp * V * ACH / 3600.
        let infiltration_conductance_w_per_k: Vec<f64> = zone_volumes_m3
            .iter()
            .map(|v| 1.2 * 1005.0 * v * infiltration_ach / 3600.0)
            .collect();

        // Inter-zone adjacency list.
        let mut neighbors: Vec<Vec<(usize, f64)>> = vec![vec![]; n];
        let mut sum_interzone_k = vec![0.0; n];
        let uid_to_idx: std::collections::HashMap<&str, usize> = zone_uids
            .iter()
            .enumerate()
            .map(|(i, uid)| (uid.as_str(), i))
            .collect();

        for e in network.interzone_conductances() {
            let Some(&i) = uid_to_idx.get(e.zone_a.as_str()) else {
                continue;
            };
            let Some(&j) = uid_to_idx.get(e.zone_b.as_str()) else {
                continue;
            };
            if i == j {
                continue;
            }
            neighbors[i].push((j, e.conductance_w_per_k));
            neighbors[j].push((i, e.conductance_w_per_k));
            sum_interzone_k[i] += e.conductance_w_per_k;
            sum_interzone_k[j] += e.conductance_w_per_k;
        }

        Self {
            zone_uids,
            zone_names,
            temperatures_c,
            thermal_capacity_j_per_k,
            exterior_conductance_w_per_k,
            infiltration_conductance_w_per_k,
            neighbors,
            sum_interzone_k,
        }
    }

    pub fn zone_names(&self) -> &[String] {
        &self.zone_names
    }

    pub fn zone_uids(&self) -> &[UID] {
        &self.zone_uids
    }

    pub fn temperatures_c(&self) -> &[f64] {
        &self.temperatures_c
    }

    /// Advances the model by one timestep.
    ///
    /// - `outdoor_temp_c`: boundary temperature for exterior + infiltration conductances
    /// - `gains_w`: per-zone gains (W), same ordering as `zone_names()`
    /// - `hvac`: ideal-loads controller (setpoints)
    /// - `dt_s`: timestep in seconds
    pub fn step(
        &mut self,
        outdoor_temp_c: f64,
        gains_w: &[f64],
        hvac: &HvacIdealLoads,
        dt_s: f64,
    ) -> Result<MultiZoneStepResult> {
        let n = self.zone_uids.len();
        anyhow::ensure!(gains_w.len() == n, "gains length mismatch");
        anyhow::ensure!(dt_s > 0.0, "dt must be positive");

        // Precompute diagonal terms and base RHS (HVAC excluded).
        let mut diag = vec![0.0; n];
        let mut rhs_base = vec![0.0; n];
        for i in 0..n {
            let c = self.thermal_capacity_j_per_k[i].max(0.0);
            let k_out =
                self.exterior_conductance_w_per_k[i] + self.infiltration_conductance_w_per_k[i];
            let c_over_dt = if c > 0.0 { c / dt_s } else { 0.0 };

            diag[i] = c_over_dt + k_out + self.sum_interzone_k[i];
            rhs_base[i] = c_over_dt * self.temperatures_c[i] + k_out * outdoor_temp_c + gains_w[i];
        }

        let mut control = vec![Control::None; n];
        let mut setpoints = vec![0.0; n];

        // Iteratively determine which zones need setpoint control, while accounting for
        // inter-zone coupling.
        for _ in 0..20 {
            let temps =
                solve_with_controls(&diag, &rhs_base, &self.neighbors, &control, &setpoints)?;

            let mut new_control = vec![Control::None; n];
            let mut new_setpoints = vec![0.0; n];

            // For each zone, compute the local "free-floating" temperature at t+dt assuming
            // HVAC is off in this zone, but neighbors remain at their solved temperatures.
            for i in 0..n {
                let c = self.thermal_capacity_j_per_k[i].max(0.0);
                let k_out =
                    self.exterior_conductance_w_per_k[i] + self.infiltration_conductance_w_per_k[i];
                let c_over_dt = if c > 0.0 { c / dt_s } else { 0.0 };
                let denom = c_over_dt + k_out + self.sum_interzone_k[i];
                let mut neighbor_term = 0.0;
                for &(j, k) in &self.neighbors[i] {
                    neighbor_term += k * temps[j];
                }
                let t_free = if denom > 1e-14 {
                    (c_over_dt * self.temperatures_c[i]
                        + k_out * outdoor_temp_c
                        + gains_w[i]
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
                // Final solve with stable control set.
                let final_temps =
                    solve_with_controls(&diag, &rhs_base, &self.neighbors, &control, &setpoints)?;

                let mut zone_heating_w = vec![0.0; n];
                let mut zone_cooling_w = vec![0.0; n];

                for i in 0..n {
                    if control[i] == Control::None {
                        continue;
                    }
                    let c = self.thermal_capacity_j_per_k[i].max(0.0);
                    let k_out = self.exterior_conductance_w_per_k[i]
                        + self.infiltration_conductance_w_per_k[i];
                    let c_over_dt = if c > 0.0 { c / dt_s } else { 0.0 };

                    let mut coupling = 0.0;
                    for &(j, k) in &self.neighbors[i] {
                        coupling += k * (final_temps[i] - final_temps[j]);
                    }

                    // Q_hvac = C/dt*(T_new - T_prev) + K_out*(T_new - T_out) + ΣKij*(T_new - Tj) - gains
                    let q_hvac = c_over_dt * (final_temps[i] - self.temperatures_c[i])
                        + k_out * (final_temps[i] - outdoor_temp_c)
                        + coupling
                        - gains_w[i];

                    if q_hvac > 0.0 {
                        zone_heating_w[i] = q_hvac;
                    } else {
                        zone_cooling_w[i] = -q_hvac;
                    }
                }

                self.temperatures_c.clone_from(&final_temps);

                return Ok(MultiZoneStepResult {
                    zone_uids: self.zone_uids.clone(),
                    zone_names: self.zone_names.clone(),
                    zone_temperatures_c: final_temps,
                    zone_heating_w,
                    zone_cooling_w,
                });
            }

            control = new_control;
            setpoints = new_setpoints;
        }

        anyhow::bail!("Multi-zone control iteration did not converge");
    }
}

fn solve_with_controls(
    diag: &[f64],
    rhs_base: &[f64],
    neighbors: &[Vec<(usize, f64)>],
    control: &[Control],
    setpoints: &[f64],
) -> Result<Vec<f64>> {
    let n = diag.len();
    let mut is_fixed = vec![false; n];
    let mut fixed_temp = vec![0.0; n];
    for i in 0..n {
        if control[i] != Control::None {
            is_fixed[i] = true;
            fixed_temp[i] = setpoints[i];
        }
    }

    let mut unknown = Vec::new();
    let mut pos = vec![usize::MAX; n];
    for i in 0..n {
        if !is_fixed[i] {
            pos[i] = unknown.len();
            unknown.push(i);
        }
    }

    let mut temps = vec![0.0; n];
    for i in 0..n {
        if is_fixed[i] {
            temps[i] = fixed_temp[i];
        }
    }

    if unknown.is_empty() {
        return Ok(temps);
    }

    let m = unknown.len();
    let mut a = vec![vec![0.0; m]; m];
    let mut b = vec![0.0; m];

    for (row_idx, &i) in unknown.iter().enumerate() {
        a[row_idx][row_idx] = diag[i];
        let mut rhs = rhs_base[i];
        for &(j, k) in &neighbors[i] {
            if is_fixed[j] {
                rhs += k * fixed_temp[j];
            } else {
                let col_idx = pos[j];
                a[row_idx][col_idx] -= k;
            }
        }
        b[row_idx] = rhs;
    }

    let x = solve_dense(a, b)?;
    for (row_idx, &i) in unknown.iter().enumerate() {
        temps[i] = x[row_idx];
    }

    Ok(temps)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sim::energy::boundary::ThermalBoundaries;
    use crate::sim::energy::config::ThermalConfig;
    use crate::sim::index::SurfaceIndex;
    use crate::{Solid, Zone};

    #[test]
    fn test_interzone_coupling_transfers_heat() {
        // Two zones sharing a face; inject gains into zone 0 and verify zone 1 warms up.
        let s0 = Solid::from_box(1.0, 1.0, 1.0, None, "s0").unwrap();
        let s1 = Solid::from_box(1.0, 1.0, 1.0, Some((1.0, 0.0, 0.0)), "s1").unwrap();
        let z0 = Zone::new("z0", vec![s0]).unwrap();
        let z1 = Zone::new("z1", vec![s1]).unwrap();
        let building = Building::new("b", vec![z0, z1]).unwrap();

        let index = SurfaceIndex::new(&building);
        let boundaries = ThermalBoundaries::classify(&building, &index);
        let mut config = ThermalConfig::new();
        config.infiltration_ach = 0.0;
        config.outdoor_temperature = 20.0;
        config.indoor_temperature = 20.0;

        let network = ThermalNetwork::build(&building, &config, &index, &boundaries);
        let mut model = MultiZoneAirModel::new(&building, &network, config.infiltration_ach, 20.0);

        // Disable HVAC by using an enormous deadband.
        let hvac = HvacIdealLoads::with_setpoints(-1e9, 1e9);
        let gains = vec![1000.0, 0.0];

        let res = model.step(20.0, &gains, &hvac, 3600.0).unwrap();
        assert!(
            res.zone_temperatures_c[1] > 20.0,
            "Zone 1 should warm via inter-zone coupling"
        );
    }

    #[test]
    fn test_no_coupling_no_temperature_change_in_other_zone() {
        // Same as above, but with a gap so there is no interface between zones.
        let s0 = Solid::from_box(1.0, 1.0, 1.0, None, "s0").unwrap();
        let s1 = Solid::from_box(1.0, 1.0, 1.0, Some((3.0, 0.0, 0.0)), "s1").unwrap();
        let z0 = Zone::new("z0", vec![s0]).unwrap();
        let z1 = Zone::new("z1", vec![s1]).unwrap();
        let building = Building::new("b", vec![z0, z1]).unwrap();

        let index = SurfaceIndex::new(&building);
        let boundaries = ThermalBoundaries::classify(&building, &index);
        let mut config = ThermalConfig::new();
        config.infiltration_ach = 0.0;
        config.outdoor_temperature = 20.0;
        config.indoor_temperature = 20.0;

        let network = ThermalNetwork::build(&building, &config, &index, &boundaries);
        let mut model = MultiZoneAirModel::new(&building, &network, config.infiltration_ach, 20.0);

        let hvac = HvacIdealLoads::with_setpoints(-1e9, 1e9);
        let gains = vec![1000.0, 0.0];

        let res = model.step(20.0, &gains, &hvac, 3600.0).unwrap();
        assert!(
            (res.zone_temperatures_c[1] - 20.0).abs() < 1e-6,
            "Without inter-zone coupling, zone 1 should remain at the boundary temperature"
        );
    }
}
