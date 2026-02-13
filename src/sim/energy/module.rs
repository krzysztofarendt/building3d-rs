use anyhow::Result;
use std::collections::HashSet;

use crate::sim::coupling::{
    InternalGainsWPerZone, InternalGainsWTotal, OutdoorAirTemperatureC,
    OutdoorWindSpeedMPerS, ShortwaveAbsorbedWPerPolygon, ShortwaveTransmittedWPerPolygon,
    ShortwaveTransmittedWPerZone,
};
use super::convection::{exterior_convection_h, interior_convection_h};
use crate::sim::framework::{Bus, SimContext, SimModule};
use crate::sim::heat_transfer::{BoundaryCondition, FvmWallSolver, build_1d_mesh};

use super::boundary::ThermalBoundaries;
use super::config::{InternalMassBoundary, ThermalConfig};
use super::hvac::HvacIdealLoads;
use super::network::{MultiZoneAirModel, MultiZoneEnvelopeRcModel, ThermalNetwork};
use super::view_factors::{
    InternalMassInfo, SurfaceHandle, ViewFactorData, compute_building_view_factors,
    compute_per_surface_mrt, linearized_h_rad,
};
use crate::sim::materials::MaterialLibrary;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EnergyModelKind {
    /// Zone air nodes only (current baseline).
    AirOnly,
    /// Zone air nodes + one aggregated 2R1C envelope node per zone.
    EnvelopeRc2R1C,
}

impl Default for EnergyModelKind {
    fn default() -> Self {
        Self::AirOnly
    }
}

/// Configuration for step-based multi-zone thermal simulation.
#[derive(Debug, Clone)]
pub struct EnergyModuleConfig {
    /// Thermal envelope and model parameters.
    pub thermal: ThermalConfig,
    /// Ideal HVAC controller.
    pub hvac: HvacIdealLoads,
    /// Timestep in seconds.
    pub dt_s: f64,
    /// If true, run in steady-state mode (equivalent to zero thermal capacity).
    pub steady_state: bool,
    /// Selects the internal thermal network model.
    pub model_kind: EnergyModelKind,
    /// Optional material library (legacy convenience wiring).
    ///
    /// If `thermal.material_library` is not set, this is copied into it during `init()`.
    /// Prefer setting `thermal.material_library` directly so U-values and capacities are
    /// resolved consistently from one source of truth.
    pub material_library: Option<MaterialLibrary>,
    /// Fallback envelope capacity per exterior area [J/(m²·K)] for surfaces that
    /// do not have `ThermalMaterial.thermal_capacity` in the material library.
    pub default_envelope_capacity_j_per_m2_k: f64,
}

impl Default for EnergyModuleConfig {
    fn default() -> Self {
        Self {
            thermal: ThermalConfig::new(),
            hvac: HvacIdealLoads::new(),
            dt_s: 3600.0,
            steady_state: false,
            model_kind: EnergyModelKind::AirOnly,
            material_library: None,
            default_envelope_capacity_j_per_m2_k: 0.0,
        }
    }
}

/// Step-based multi-zone thermal simulation module.
///
/// Inputs (via [`Bus`]):
/// - [`OutdoorAirTemperatureC`] (optional; falls back to `thermal.outdoor_temperature`)
/// - [`InternalGainsWPerZone`] or [`InternalGainsWTotal`] (optional; default 0)
/// - [`ShortwaveTransmittedWPerZone`] (optional; treated as solar shortwave gains)
///
/// Outputs (via [`Bus`]):
/// - [`MultiZoneStepResult`] for the latest step
pub struct EnergyModule {
    config: EnergyModuleConfig,
    model: Option<EnergyModel>,
    step_index: usize,
    zone_volumes_m3: Vec<f64>,
    total_volume_m3: f64,
    zone_uids: Vec<crate::UID>,
    polygon_uid_to_zone_idx: std::collections::HashMap<crate::UID, usize>,
    boundaries: Option<ThermalBoundaries>,
    fvm_walls: Vec<FvmExteriorWall>,
    fvm_polygon_uids: HashSet<crate::UID>,
    ground_ua_by_zone_w_per_k: Vec<f64>,
    internal_mass_surfaces: Vec<FvmInternalMassSurface>,
    steady_state_exterior_surfaces: Vec<SteadyStateExteriorSurface>,
    view_factor_data: Option<ViewFactorData>,
}

enum EnergyModel {
    Air(MultiZoneAirModel),
    EnvelopeRc(MultiZoneEnvelopeRcModel),
}

struct FvmExteriorWall {
    polygon_uid: crate::UID,
    zone_idx: usize,
    area_m2: f64,
    is_ground_coupled: bool,
    h_in_w_per_m2_k: f64,
    h_out_w_per_m2_k: f64,
    cos_tilt: f64,
    h_min_iso_interior: f64,
    solver: FvmWallSolver,
    /// Cached interior surface temperature from the last solver step.
    cached_interior_surface_temp_c: f64,
}

struct FvmInternalMassSurface {
    zone_idx: usize,
    /// Index of this mass surface (for `SurfaceHandle::InternalMass`).
    mass_index: usize,
    area_m2: f64,
    boundary: InternalMassBoundary,
    h_w_per_m2_k: f64,
    cos_tilt: f64,
    h_min_iso: f64,
    solver: FvmWallSolver,
    /// Cached interior surface temperature from the last solver step.
    cached_interior_surface_temp_c: f64,
    /// Cached exterior surface temperature (only meaningful for TwoSided).
    cached_exterior_surface_temp_c: f64,
}

/// A non-FVM exterior surface (e.g. glazing) whose interior surface temperature
/// is estimated via steady-state heat balance for inclusion in MRT.
struct SteadyStateExteriorSurface {
    polygon_uid: crate::UID,
    zone_idx: usize,
    area_m2: f64,
    u_value_w_per_m2_k: f64,
    h_in_w_per_m2_k: f64,
    cos_tilt: f64,
}

impl EnergyModule {
    pub fn new(config: EnergyModuleConfig) -> Self {
        Self {
            config,
            model: None,
            step_index: 0,
            zone_volumes_m3: vec![],
            total_volume_m3: 0.0,
            zone_uids: vec![],
            polygon_uid_to_zone_idx: std::collections::HashMap::new(),
            boundaries: None,
            fvm_walls: vec![],
            fvm_polygon_uids: HashSet::new(),
            ground_ua_by_zone_w_per_k: vec![],
            internal_mass_surfaces: vec![],
            steady_state_exterior_surfaces: vec![],
            view_factor_data: None,
        }
    }

    fn gains_by_zone_split(&self, bus: &Bus) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
        let n = self.zone_uids.len();
        let mut air_gains = vec![0.0; n];
        let mut env_gains = vec![0.0; n];
        let mut internal_mass_sources = vec![0.0; n];

        let internal_by_zone = bus
            .get::<InternalGainsWPerZone>()
            .map(|g| &g.watts_by_zone_uid);

        let internal_total = bus.get::<InternalGainsWTotal>().map(|g| g.0).unwrap_or(0.0);

        let solar_by_zone = bus
            .get::<ShortwaveTransmittedWPerZone>()
            .map(|g| &g.watts_by_zone_uid);

        let use_internal_mass_sources = self.config.thermal.use_surface_aware_solar_distribution
            && !self.internal_mass_surfaces.is_empty();
        let f_internal_mass = self
            .config
            .thermal
            .internal_gains_to_mass_fraction
            .clamp(0.0, 1.0);
        let f_solar_air = self
            .config
            .thermal
            .transmitted_solar_to_air_fraction
            .clamp(0.0, 1.0);

        // Apply internal + transmitted shortwave to the air node, optionally splitting a
        // portion to internal mass slabs for lagged release.
        for (i, air_gain) in air_gains.iter_mut().enumerate() {
            let uid = &self.zone_uids[i];

            let internal_i = if let Some(map) = internal_by_zone {
                map.get(uid).cloned().unwrap_or(0.0)
            } else if self.total_volume_m3 > 1e-14 {
                internal_total * (self.zone_volumes_m3[i] / self.total_volume_m3)
            } else {
                0.0
            };

            let solar_i = solar_by_zone
                .and_then(|map| map.get(uid).cloned())
                .unwrap_or(0.0);

            if use_internal_mass_sources {
                *air_gain += internal_i * (1.0 - f_internal_mass) + solar_i * f_solar_air;
                internal_mass_sources[i] +=
                    internal_i * f_internal_mass + solar_i * (1.0 - f_solar_air);
            } else {
                *air_gain += internal_i + solar_i;
            }
        }

        // When SolarInteriorDistributionModule has published per-polygon
        // transmitted solar, use that (aggregated per zone) for internal mass sources.
        // This replaces the bulk per-zone transmitted solar for the mass path,
        // giving the distribution module control over beam vs diffuse allocation.
        if use_internal_mass_sources
            && let Some(per_poly) = bus.get::<ShortwaveTransmittedWPerPolygon>()
        {
            // Re-aggregate per zone from per-polygon data.
            let mut per_zone_from_poly = vec![0.0_f64; n];
            for (polygon_uid, w) in &per_poly.watts_by_polygon_uid {
                if *w == 0.0 {
                    continue;
                }
                if let Some(&zone_idx) = self.polygon_uid_to_zone_idx.get(polygon_uid) {
                    per_zone_from_poly[zone_idx] += *w;
                }
            }
            // Replace internal_mass_sources solar component with the distributed values.
            for (i, ims) in internal_mass_sources.iter_mut().enumerate() {
                let uid = &self.zone_uids[i];
                let solar_zone_orig = solar_by_zone
                    .and_then(|map| map.get(uid).cloned())
                    .unwrap_or(0.0);
                // Remove original solar contribution and add the distributed one.
                let solar_to_mass_orig = solar_zone_orig * (1.0 - f_solar_air);
                let solar_to_mass_new = per_zone_from_poly[i] * (1.0 - f_solar_air);
                *ims += solar_to_mass_new - solar_to_mass_orig;
            }
        }

        // Per-polygon absorbed shortwave: split between envelope (exterior) and air (non-exterior).
        let is_exterior = |uid: &crate::UID| {
            self.boundaries
                .as_ref()
                .map(|b| b.is_exterior(uid))
                .unwrap_or(false)
        };

        if let Some(abs) = bus.get::<ShortwaveAbsorbedWPerPolygon>() {
            for (polygon_uid, w) in &abs.watts_by_polygon_uid {
                if *w == 0.0 {
                    continue;
                }
                // If this polygon is modeled via an FVM wall solver, its absorbed shortwave
                // is applied as a boundary flux to the wall (not as an instantaneous gain).
                if self.fvm_polygon_uids.contains(polygon_uid) {
                    continue;
                }
                if let Some(&zone_idx) = self.polygon_uid_to_zone_idx.get(polygon_uid) {
                    if is_exterior(polygon_uid) {
                        env_gains[zone_idx] += *w;
                    } else {
                        air_gains[zone_idx] += *w;
                    }
                }
            }
        }

        (air_gains, env_gains, internal_mass_sources)
    }
}

impl SimModule for EnergyModule {
    fn name(&self) -> &'static str {
        "energy"
    }

    fn init(&mut self, ctx: &SimContext, _bus: &mut Bus) -> Result<()> {
        // Backward compatible wiring: allow the module-level material library field to
        // populate `ThermalConfig.material_library` when the latter is unset.
        if self.config.thermal.material_library.is_none() {
            self.config.thermal.material_library = self.config.material_library.clone();
        }

        let boundaries = ThermalBoundaries::classify(ctx.building, ctx.surface_index);
        self.boundaries = Some(boundaries.clone());

        let zones = ctx.building.zones();
        self.zone_volumes_m3 = zones.iter().map(|z| z.volume()).collect();
        self.total_volume_m3 = self.zone_volumes_m3.iter().sum();
        self.zone_uids = zones.iter().map(|z| z.uid.clone()).collect();

        self.polygon_uid_to_zone_idx.clear();
        let uid_to_idx: std::collections::HashMap<&str, usize> = self
            .zone_uids
            .iter()
            .enumerate()
            .map(|(i, uid)| (uid.as_str(), i))
            .collect();
        for s in &ctx.surface_index.surfaces {
            if let Some(&idx) = uid_to_idx.get(s.zone_uid.as_str()) {
                self.polygon_uid_to_zone_idx
                    .insert(s.polygon_uid.clone(), idx);
            }
        }

        self.fvm_walls.clear();
        self.fvm_polygon_uids.clear();
        if self.config.thermal.use_fvm_walls {
            for s in &ctx.surface_index.surfaces {
                if !boundaries.is_exterior(&s.polygon_uid) {
                    continue;
                }
                if s.area_m2 <= 0.0 {
                    continue;
                }
                if looks_like_glazing(&s.path, &self.config.thermal) {
                    continue;
                }
                if self
                    .config
                    .thermal
                    .has_u_value_override_for_surface(&s.polygon_uid, &s.path)
                {
                    continue;
                }
                if self.config.thermal.resolve_construction(&s.path).is_none() {
                    continue;
                }

                let is_ground_coupled = if self.config.thermal.ground_temperature_c.is_some()
                    && self
                        .config
                        .thermal
                        .ground_surface_patterns
                        .iter()
                        .any(|p| s.path.contains(p.as_str()))
                {
                    ctx.building
                        .get_polygon(&s.path)
                        .is_some_and(|poly| poly.vn.dz <= -0.5)
                } else {
                    false
                };
                if is_ground_coupled {
                    continue;
                }

                let Some(&zone_idx) = self.polygon_uid_to_zone_idx.get(&s.polygon_uid) else {
                    continue;
                };
                let Some(construction) = self.config.thermal.resolve_construction(&s.path) else {
                    continue;
                };

                let h_si_iso = if construction.r_si > 0.0 {
                    1.0 / construction.r_si
                } else {
                    1.0 / 0.13
                }
                .max(1e-9);
                let h_in = if self.config.thermal.interior_heat_transfer_coeff_w_per_m2_k > 0.0 {
                    self.config
                        .thermal
                        .interior_heat_transfer_coeff_w_per_m2_k
                        .max(h_si_iso)
                } else {
                    h_si_iso
                };
                let h_out = if construction.r_se > 0.0 {
                    1.0 / construction.r_se
                } else {
                    1.0 / 0.04
                }
                .max(1e-9);

                let cos_tilt = ctx
                    .building
                    .get_polygon(&s.path)
                    .map(|p| p.vn.dz)
                    .unwrap_or(0.0);

                let mesh = build_1d_mesh(construction, s.area_m2);
                let solver = FvmWallSolver::new(mesh, self.config.thermal.indoor_temperature);
                let init_temp = solver.interior_surface_temp();
                self.fvm_walls.push(FvmExteriorWall {
                    polygon_uid: s.polygon_uid.clone(),
                    zone_idx,
                    area_m2: s.area_m2,
                    is_ground_coupled: false,
                    h_in_w_per_m2_k: h_in,
                    h_out_w_per_m2_k: h_out,
                    cos_tilt,
                    h_min_iso_interior: h_si_iso,
                    solver,
                    cached_interior_surface_temp_c: init_temp,
                });
                self.fvm_polygon_uids.insert(s.polygon_uid.clone());
            }
        }

        self.ground_ua_by_zone_w_per_k = vec![0.0; self.zone_uids.len()];
        if self.config.thermal.ground_temperature_c.is_some() {
            for s in &ctx.surface_index.surfaces {
                if !boundaries.is_exterior(&s.polygon_uid) {
                    continue;
                }
                if self.fvm_polygon_uids.contains(&s.polygon_uid) {
                    continue;
                }
                if !self
                    .config
                    .thermal
                    .ground_surface_patterns
                    .iter()
                    .any(|p| s.path.contains(p.as_str()))
                {
                    continue;
                }
                let Some(poly) = ctx.building.get_polygon(&s.path) else {
                    continue;
                };
                if poly.vn.dz > -0.5 {
                    continue;
                }

                let Some(&zone_idx) = self.polygon_uid_to_zone_idx.get(&s.polygon_uid) else {
                    continue;
                };
                let u = self
                    .config
                    .thermal
                    .resolve_u_value_for_surface(&s.polygon_uid, &s.path);
                self.ground_ua_by_zone_w_per_k[zone_idx] += u * s.area_m2;
            }
        }

        // Collect non-FVM exterior surfaces (e.g. glazing, U-value-only surfaces) for
        // steady-state interior surface temperature in MRT.
        self.steady_state_exterior_surfaces.clear();
        for s in &ctx.surface_index.surfaces {
            if !boundaries.is_exterior(&s.polygon_uid) {
                continue;
            }
            if s.area_m2 <= 0.0 {
                continue;
            }
            // Skip surfaces already modeled by FVM solvers.
            if self.fvm_polygon_uids.contains(&s.polygon_uid) {
                continue;
            }
            // Skip ground-coupled surfaces (they exchange with ground, not outdoor air).
            if self.config.thermal.ground_temperature_c.is_some()
                && self
                    .config
                    .thermal
                    .ground_surface_patterns
                    .iter()
                    .any(|p| s.path.contains(p.as_str()))
                && let Some(poly) = ctx.building.get_polygon(&s.path)
                && poly.vn.dz <= -0.5
            {
                continue;
            }
            let Some(&zone_idx) = self.polygon_uid_to_zone_idx.get(&s.polygon_uid) else {
                continue;
            };
            let u = self
                .config
                .thermal
                .resolve_u_value_for_surface(&s.polygon_uid, &s.path);
            if u <= 0.0 {
                continue;
            }
            let h_si_iso = 1.0 / 0.13_f64;
            let h_in = if self.config.thermal.interior_heat_transfer_coeff_w_per_m2_k > 0.0 {
                self.config
                    .thermal
                    .interior_heat_transfer_coeff_w_per_m2_k
                    .max(h_si_iso)
            } else {
                h_si_iso
            };
            let cos_tilt = ctx
                .building
                .get_polygon(&s.path)
                .map(|p| p.vn.dz)
                .unwrap_or(0.0);
            self.steady_state_exterior_surfaces
                .push(SteadyStateExteriorSurface {
                    polygon_uid: s.polygon_uid.clone(),
                    zone_idx,
                    area_m2: s.area_m2,
                    u_value_w_per_m2_k: u,
                    h_in_w_per_m2_k: h_in,
                    cos_tilt,
                });
        }

        self.internal_mass_surfaces.clear();
        if !self.config.thermal.internal_mass_surfaces.is_empty() {
            let mut mass_index_counter = 0_usize;
            for (zone_idx, zone) in ctx.building.zones().iter().enumerate() {
                for m in &self.config.thermal.internal_mass_surfaces {
                    if m.area_m2 <= 0.0 {
                        continue;
                    }
                    if !zone.name.contains(m.zone_path_pattern.as_str()) {
                        continue;
                    }
                    let h_si_iso = if m.construction.r_si > 0.0 {
                        1.0 / m.construction.r_si
                    } else {
                        3.0
                    }
                    .max(1e-9);
                    let h = if self.config.thermal.interior_heat_transfer_coeff_w_per_m2_k > 0.0 {
                        self.config
                            .thermal
                            .interior_heat_transfer_coeff_w_per_m2_k
                            .max(h_si_iso)
                    } else {
                        h_si_iso
                    };
                    let this_mass_index = mass_index_counter;
                    mass_index_counter += 1;
                    let mesh = build_1d_mesh(&m.construction, m.area_m2);
                    let solver = FvmWallSolver::new(mesh, self.config.thermal.indoor_temperature);
                    let init_temp_in = solver.interior_surface_temp();
                    let init_temp_out = solver.exterior_surface_temp();
                    self.internal_mass_surfaces.push(FvmInternalMassSurface {
                        zone_idx,
                        mass_index: this_mass_index,
                        area_m2: m.area_m2,
                        boundary: m.boundary,
                        h_w_per_m2_k: h,
                        cos_tilt: m.cos_tilt,
                        h_min_iso: h_si_iso,
                        solver,
                        cached_interior_surface_temp_c: init_temp_in,
                        cached_exterior_surface_temp_c: init_temp_out,
                    });
                }
            }
        }

        let network = if self.config.thermal.use_fvm_walls && !self.fvm_polygon_uids.is_empty() {
            ThermalNetwork::build_with_ignored_exterior_polygons(
                ctx.building,
                &self.config.thermal,
                ctx.surface_index,
                &boundaries,
                &self.fvm_polygon_uids,
            )
        } else {
            ThermalNetwork::build(
                ctx.building,
                &self.config.thermal,
                ctx.surface_index,
                &boundaries,
            )
        };

        let mut cap_j_per_m3_k = if self.config.steady_state {
            0.0
        } else {
            self.config.thermal.thermal_capacity_j_per_m3_k
        };

        // When FVM walls and/or internal mass slabs are present, all structural
        // mass is explicitly modeled. Use pure air capacity (rho * cp) instead of
        // the lumped volumetric capacity to avoid double-counting.
        if ((!self.fvm_walls.is_empty()) || (!self.internal_mass_surfaces.is_empty()))
            && self.total_volume_m3 > 0.0
            && cap_j_per_m3_k > 0.0
        {
            cap_j_per_m3_k = 1.2 * 1005.0; // pure air: rho*cp [J/(m³·K)]
        }

        self.model = Some(match self.config.model_kind {
            EnergyModelKind::AirOnly => EnergyModel::Air(MultiZoneAirModel::new(
                ctx.building,
                &network,
                self.config.thermal.infiltration_ach,
                cap_j_per_m3_k,
                self.config.thermal.indoor_temperature,
            )),
            EnergyModelKind::EnvelopeRc2R1C => {
                let default_env_cap = if self.config.steady_state {
                    0.0
                } else {
                    self.config.default_envelope_capacity_j_per_m2_k
                };

                EnergyModel::EnvelopeRc(MultiZoneEnvelopeRcModel::new_with_thermal_config(
                    ctx.building,
                    &network,
                    ctx.surface_index,
                    &boundaries,
                    &self.config.thermal,
                    self.config.thermal.infiltration_ach,
                    cap_j_per_m3_k,
                    default_env_cap,
                    self.config.thermal.indoor_temperature,
                ))
            }
        });

        // Compute per-zone view factors when enabled.
        self.view_factor_data = if self.config.thermal.use_view_factor_radiation {
            let mass_infos: Vec<InternalMassInfo> = self
                .internal_mass_surfaces
                .iter()
                .map(|m| InternalMassInfo {
                    index: m.mass_index,
                    zone_uid: self.zone_uids[m.zone_idx].clone(),
                    cos_tilt: m.cos_tilt,
                    area_m2: m.area_m2,
                })
                .collect();
            Some(compute_building_view_factors(
                ctx.building,
                ctx.surface_index,
                &boundaries,
                &mass_infos,
                self.config.thermal.view_factor_rays_per_surface,
            ))
        } else {
            None
        };

        Ok(())
    }

    fn step(&mut self, _ctx: &SimContext, bus: &mut Bus) -> Result<()> {
        let outdoor_temp_c = bus
            .get::<OutdoorAirTemperatureC>()
            .map(|t| t.0)
            .unwrap_or(self.config.thermal.outdoor_temperature);
        let wind_speed = bus
            .get::<OutdoorWindSpeedMPerS>()
            .map(|w| w.0)
            .unwrap_or(0.0);

        let (mut air_gains, mut env_gains, mut internal_mass_sources) = self.gains_by_zone_split(bus);
        if let Some(tg) = self.config.thermal.ground_temperature_c {
            let dt = tg - outdoor_temp_c;
            for (i, ua) in self.ground_ua_by_zone_w_per_k.iter().enumerate() {
                if let Some(g) = env_gains.get_mut(i) {
                    *g += ua * dt;
                }
            }
        }

        let Some(model) = self.model.as_mut() else {
            anyhow::bail!("EnergyModule not initialized");
        };

        let air_temps: Vec<f64> = match model {
            EnergyModel::Air(m) => m.temperatures_c().to_vec(),
            EnergyModel::EnvelopeRc(m) => m.air_temperatures_c().to_vec(),
        };
        // ── Per-surface MRT (view factors) or area-weighted MRT (legacy) ──
        let (vf_mrt_map, vf_h_r) = if let Some(vf) = &self.view_factor_data {
            let mut surf_temps: std::collections::HashMap<SurfaceHandle, f64> =
                std::collections::HashMap::new();
            for w in &self.fvm_walls {
                surf_temps.insert(
                    SurfaceHandle::Polygon(w.polygon_uid.clone()),
                    w.cached_interior_surface_temp_c,
                );
            }
            for m in &self.internal_mass_surfaces {
                surf_temps.insert(
                    SurfaceHandle::InternalMass {
                        index: m.mass_index,
                    },
                    m.cached_interior_surface_temp_c,
                );
            }
            for ss in &self.steady_state_exterior_surfaces {
                let t_air_zone = air_temps
                    .get(ss.zone_idx)
                    .copied()
                    .unwrap_or(self.config.thermal.indoor_temperature);
                let ratio =
                    ss.u_value_w_per_m2_k / (ss.u_value_w_per_m2_k + ss.h_in_w_per_m2_k);
                let t_surf = t_air_zone - ratio * (t_air_zone - outdoor_temp_c);
                surf_temps.insert(SurfaceHandle::Polygon(ss.polygon_uid.clone()), t_surf);
            }
            let t_mean = air_temps.first().copied().unwrap_or(20.0);
            let mrt_map = compute_per_surface_mrt(vf, &surf_temps);
            let h_r = linearized_h_rad(self.config.thermal.interior_emissivity, t_mean);
            (Some(mrt_map), h_r)
        } else {
            (None, 0.0)
        };

        let mut rad_temps = air_temps.clone();
        if self.view_factor_data.is_none() && self.config.thermal.use_interior_radiative_exchange {
            let f_rad = self
                .config
                .thermal
                .interior_radiation_fraction
                .clamp(0.0, 1.0);
            let n = self.zone_uids.len();
            let mut num = vec![0.0_f64; n];
            let mut den = vec![0.0_f64; n];

            for w in &self.fvm_walls {
                if w.area_m2 <= 0.0 {
                    continue;
                }
                let h_rad = w.h_in_w_per_m2_k.max(1e-9) * f_rad;
                if h_rad <= 0.0 {
                    continue;
                }
                num[w.zone_idx] += h_rad * w.area_m2 * w.cached_interior_surface_temp_c;
                den[w.zone_idx] += h_rad * w.area_m2;
            }
            for m in &self.internal_mass_surfaces {
                if m.area_m2 <= 0.0 {
                    continue;
                }
                let h_rad = m.h_w_per_m2_k.max(1e-9) * f_rad;
                if h_rad <= 0.0 {
                    continue;
                }
                num[m.zone_idx] += h_rad * m.area_m2 * m.cached_interior_surface_temp_c;
                den[m.zone_idx] += h_rad * m.area_m2;
                if matches!(m.boundary, InternalMassBoundary::TwoSided) {
                    num[m.zone_idx] += h_rad * m.area_m2 * m.cached_exterior_surface_temp_c;
                    den[m.zone_idx] += h_rad * m.area_m2;
                }
            }

            // Steady-state exterior surfaces (e.g. windows): estimate interior
            // surface temp from T_surf = T_air - U/(U+h_in) * (T_air - T_out).
            for ss in &self.steady_state_exterior_surfaces {
                if ss.area_m2 <= 0.0 {
                    continue;
                }
                let t_air_zone = air_temps
                    .get(ss.zone_idx)
                    .copied()
                    .unwrap_or(self.config.thermal.indoor_temperature);
                let ratio_est =
                    ss.u_value_w_per_m2_k / (ss.u_value_w_per_m2_k + ss.h_in_w_per_m2_k);
                let t_surf_est = t_air_zone - ratio_est * (t_air_zone - outdoor_temp_c);
                let h_in_dyn = interior_convection_h(
                    &self.config.thermal.interior_convection_model,
                    t_surf_est - t_air_zone,
                    ss.cos_tilt,
                    ss.h_in_w_per_m2_k,
                )
                .max(1e-9);
                let h_rad = h_in_dyn * f_rad;
                if h_rad <= 0.0 {
                    continue;
                }
                let ratio = ss.u_value_w_per_m2_k / (ss.u_value_w_per_m2_k + h_in_dyn);
                let t_surf = t_air_zone - ratio * (t_air_zone - outdoor_temp_c);
                num[ss.zone_idx] += h_rad * ss.area_m2 * t_surf;
                den[ss.zone_idx] += h_rad * ss.area_m2;
            }

            for i in 0..n {
                if den[i] > 0.0 {
                    rad_temps[i] = num[i] / den[i];
                }
            }
        }

        if !self.internal_mass_surfaces.is_empty() {
            let n = self.zone_uids.len();
            let mut total_area_by_zone_m2 = vec![0.0; n];
            for s in &self.internal_mass_surfaces {
                if s.area_m2 <= 0.0 {
                    continue;
                }
                let w_total = internal_mass_sources
                    .get(s.zone_idx)
                    .copied()
                    .unwrap_or(0.0);
                if w_total == 0.0 {
                    continue;
                }
                if let Some(a) = total_area_by_zone_m2.get_mut(s.zone_idx) {
                    *a += s.area_m2;
                }
            }

            // When fvm_wall_solar_to_air is enabled, include FVM wall interior area
            // in the denominator (diluting mass flux) and redirect the wall share to air.
            if self.config.thermal.fvm_wall_solar_to_air
                && self.config.thermal.distribute_transmitted_solar_to_fvm_walls
            {
                let mut wall_area_by_zone = vec![0.0_f64; n];
                for w in &self.fvm_walls {
                    if w.area_m2 <= 0.0 {
                        continue;
                    }
                    if let Some(a) = wall_area_by_zone.get_mut(w.zone_idx) {
                        *a += w.area_m2;
                    }
                }
                for i in 0..n {
                    let a_walls = wall_area_by_zone[i];
                    if a_walls <= 0.0 {
                        continue;
                    }
                    let a_mass = total_area_by_zone_m2[i];
                    let a_total = a_mass + a_walls;
                    if a_total <= 0.0 {
                        continue;
                    }
                    let src = internal_mass_sources[i];
                    if src == 0.0 {
                        continue;
                    }
                    let wall_share = src * (a_walls / a_total);
                    internal_mass_sources[i] -= wall_share;
                    air_gains[i] += wall_share;
                }
            }

            for s in &mut self.internal_mass_surfaces {
                let t_air = air_temps
                    .get(s.zone_idx)
                    .copied()
                    .unwrap_or(self.config.thermal.indoor_temperature);
                // Dynamic interior h for internal mass surfaces.
                let t_surf_prev = s.solver.interior_surface_temp();
                let h_total = interior_convection_h(
                    &self.config.thermal.interior_convection_model,
                    t_surf_prev - t_air,
                    s.cos_tilt,
                    s.h_min_iso,
                )
                .max(1e-9);
                let (h_conv, h_in_total, t_eff) = if let Some(ref mrts) = vf_mrt_map {
                    // View-factor path: TARP convection + uniform linearized radiation
                    let h_conv = h_total; // h_total here is convection-only from TARP
                    let eps = self.config.thermal.interior_emissivity;
                    let h_rad = eps * vf_h_r;
                    let handle = SurfaceHandle::InternalMass { index: s.mass_index };
                    let t_mrt = mrts.get(&handle).copied().unwrap_or(t_air);
                    let h_t = h_conv + h_rad;
                    let t_eff = if h_t > 0.0 {
                        (h_conv * t_air + h_rad * t_mrt) / h_t
                    } else {
                        t_air
                    };
                    (h_conv, h_t, t_eff)
                } else if self.config.thermal.use_interior_radiative_exchange {
                    let f_rad = self
                        .config
                        .thermal
                        .interior_radiation_fraction
                        .clamp(0.0, 1.0);
                    let h_rad = h_total * f_rad;
                    let h_conv = (h_total - h_rad).max(0.0);
                    let t_rad = rad_temps
                        .get(s.zone_idx)
                        .copied()
                        .unwrap_or(self.config.thermal.indoor_temperature);
                    let t_eff = if h_total > 0.0 {
                        (h_conv * t_air + h_rad * t_rad) / h_total
                    } else {
                        t_air
                    };
                    (h_conv, h_total, t_eff)
                } else {
                    (h_total, h_total, t_air)
                };

                let w_total = internal_mass_sources
                    .get(s.zone_idx)
                    .copied()
                    .unwrap_or(0.0);
                let a_total = total_area_by_zone_m2
                    .get(s.zone_idx)
                    .copied()
                    .unwrap_or(0.0);
                let source_flux_w_per_m2 = if w_total != 0.0 && a_total > 0.0 {
                    w_total / a_total
                } else {
                    0.0
                };

                let (bc_exterior, bc_interior, direct_to_air_w) = if source_flux_w_per_m2 != 0.0 {
                    let bc_in = BoundaryCondition::ConvectiveWithFluxToDomain {
                        h: h_in_total,
                        t_fluid: t_eff,
                        heat_flux: source_flux_w_per_m2,
                    };
                    let bc_out = match s.boundary {
                        InternalMassBoundary::TwoSided => BoundaryCondition::Convective {
                            h: h_in_total,
                            t_fluid: t_eff,
                        },
                        InternalMassBoundary::OneSidedAdiabatic => {
                            BoundaryCondition::Neumann { heat_flux: 0.0 }
                        }
                    };
                    (bc_out, bc_in, 0.0)
                } else {
                    let bc_in = BoundaryCondition::Convective {
                        h: h_in_total,
                        t_fluid: t_eff,
                    };
                    let bc_out = match s.boundary {
                        InternalMassBoundary::TwoSided => BoundaryCondition::Convective {
                            h: h_in_total,
                            t_fluid: t_eff,
                        },
                        InternalMassBoundary::OneSidedAdiabatic => {
                            BoundaryCondition::Neumann { heat_flux: 0.0 }
                        }
                    };
                    (bc_out, bc_in, 0.0)
                };

                s.solver
                    .step(self.config.dt_s, &bc_exterior, &bc_interior, &[]);

                // Cache proper surface temperatures for MRT computation next substep.
                s.cached_interior_surface_temp_c =
                    s.solver.interior_surface_temperature(&bc_interior);
                if matches!(s.boundary, InternalMassBoundary::TwoSided) {
                    s.cached_exterior_surface_temp_c =
                        s.solver.exterior_surface_temperature(&bc_exterior);
                }

                let mut q_to_air_w = 0.0;
                let t_surf_in = s.cached_interior_surface_temp_c;
                if vf_mrt_map.is_some() {
                    q_to_air_w += h_in_total * (t_surf_in - t_eff) * s.area_m2;
                } else {
                    q_to_air_w += h_conv * (t_surf_in - t_air) * s.area_m2;
                }
                if matches!(s.boundary, InternalMassBoundary::TwoSided) {
                    let t_surf_out = s.cached_exterior_surface_temp_c;
                    if vf_mrt_map.is_some() {
                        q_to_air_w += h_in_total * (t_surf_out - t_eff) * s.area_m2;
                    } else {
                        q_to_air_w += h_conv * (t_surf_out - t_air) * s.area_m2;
                    }
                }
                q_to_air_w += direct_to_air_w;
                if q_to_air_w != 0.0
                    && let Some(g) = air_gains.get_mut(s.zone_idx)
                {
                    *g += q_to_air_w;
                }
            }
        }

        if !self.fvm_walls.is_empty() {
            let absorbed = bus
                .get::<ShortwaveAbsorbedWPerPolygon>()
                .map(|a| &a.watts_by_polygon_uid);

            for w in &mut self.fvm_walls {
                let t_air = air_temps
                    .get(w.zone_idx)
                    .copied()
                    .unwrap_or(self.config.thermal.indoor_temperature);
                // Dynamic interior h: use current interior surface temp for dT.
                let t_surf_prev = w.solver.interior_surface_temp();
                let h_in = interior_convection_h(
                    &self.config.thermal.interior_convection_model,
                    t_surf_prev - t_air,
                    w.cos_tilt,
                    w.h_min_iso_interior,
                )
                .max(1e-9);
                // Dynamic exterior h: use current exterior surface temp for dT.
                let t_surf_ext_prev = w.solver.exterior_surface_temp();
                let h_out = exterior_convection_h(
                    &self.config.thermal.exterior_convection_model,
                    w.h_out_w_per_m2_k,
                    t_surf_ext_prev - outdoor_temp_c,
                    w.cos_tilt,
                    wind_speed,
                )
                .max(1e-9);
                let (h_in_conv, h_in_total, t_eff) = if let Some(ref mrts) = vf_mrt_map {
                    // View-factor path: TARP convection + uniform linearized radiation
                    let h_conv = h_in; // h_in is convection-only from TARP
                    let eps = self.config.thermal.interior_emissivity;
                    let h_rad = eps * vf_h_r;
                    let handle = SurfaceHandle::Polygon(w.polygon_uid.clone());
                    let t_mrt = mrts.get(&handle).copied().unwrap_or(t_air);
                    let h_t = h_conv + h_rad;
                    let t_eff = if h_t > 0.0 {
                        (h_conv * t_air + h_rad * t_mrt) / h_t
                    } else {
                        t_air
                    };
                    (h_conv, h_t, t_eff)
                } else if self.config.thermal.use_interior_radiative_exchange {
                    let f_rad = self
                        .config
                        .thermal
                        .interior_radiation_fraction
                        .clamp(0.0, 1.0);
                    let h_rad = h_in * f_rad;
                    let h_conv = (h_in - h_rad).max(0.0);
                    let t_rad = rad_temps
                        .get(w.zone_idx)
                        .copied()
                        .unwrap_or(self.config.thermal.indoor_temperature);
                    let t_eff = if h_in > 0.0 {
                        (h_conv * t_air + h_rad * t_rad) / h_in
                    } else {
                        t_air
                    };
                    (h_conv, h_in, t_eff)
                } else {
                    (h_in, h_in, t_air)
                };
                let absorbed_w = absorbed
                    .and_then(|m| m.get(&w.polygon_uid).copied())
                    .unwrap_or(0.0);
                let heat_flux_w_per_m2 = if w.area_m2 > 0.0 {
                    absorbed_w / w.area_m2
                } else {
                    0.0
                };

                let bc_in = BoundaryCondition::Convective {
                    h: h_in_total,
                    t_fluid: t_eff,
                };
                let bc_out = if w.is_ground_coupled {
                    BoundaryCondition::Convective {
                        h: h_out,
                        t_fluid: self
                            .config
                            .thermal
                            .ground_temperature_c
                            .unwrap_or(outdoor_temp_c),
                    }
                } else if heat_flux_w_per_m2 != 0.0 {
                    BoundaryCondition::ConvectiveWithFlux {
                        h: h_out,
                        t_fluid: outdoor_temp_c,
                        heat_flux: heat_flux_w_per_m2,
                    }
                } else {
                    BoundaryCondition::Convective {
                        h: h_out,
                        t_fluid: outdoor_temp_c,
                    }
                };

                w.solver.step(self.config.dt_s, &bc_out, &bc_in, &[]);
                w.cached_interior_surface_temp_c =
                    w.solver.interior_surface_temperature(&bc_in);
                let t_surf = w.cached_interior_surface_temp_c;
                let q_in_w = if vf_mrt_map.is_some() {
                    h_in_total * (t_surf - t_eff) * w.area_m2
                } else {
                    h_in_conv * (t_surf - t_air) * w.area_m2
                };
                if let Some(g) = air_gains.get_mut(w.zone_idx) {
                    *g += q_in_w;
                }
            }
        }

        let result = match model {
            EnergyModel::Air(m) => {
                let mut gains = air_gains;
                for (i, v) in env_gains.into_iter().enumerate() {
                    gains[i] += v;
                }
                m.step(outdoor_temp_c, &gains, &self.config.hvac, self.config.dt_s)?
            }
            EnergyModel::EnvelopeRc(m) => m.step(
                outdoor_temp_c,
                &air_gains,
                &env_gains,
                &self.config.hvac,
                self.config.dt_s,
            )?,
        };
        bus.put(result.clone());

        self.step_index += 1;
        Ok(())
    }
}

fn looks_like_glazing(path: &str, thermal: &ThermalConfig) -> bool {
    if let Some(lib) = thermal.material_library.as_ref()
        && let Some(mat) = lib.lookup(path)
        && mat.is_glazing
    {
        return true;
    }
    let p = path.to_ascii_lowercase();
    p.contains("window") || p.contains("glazing") || p.contains("glass")
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sim::energy::network::MultiZoneStepResult;
    use crate::sim::framework::SimContext;
    use crate::sim::index::SurfaceIndex;
    use crate::{Building, Solid, Zone};

    #[test]
    fn test_energy_module_name() {
        let module = EnergyModule::new(EnergyModuleConfig::default());
        assert_eq!(module.name(), "energy");
    }

    #[test]
    fn test_energy_module_config_default_trait() {
        let cfg: EnergyModuleConfig = Default::default();
        assert_eq!(cfg.model_kind, EnergyModelKind::AirOnly);
        assert!((cfg.dt_s - 3600.0).abs() < 1e-12);
    }

    #[test]
    fn test_energy_module_step_without_init_errors() {
        let s = Solid::from_box(1.0, 1.0, 1.0, None, "s").unwrap();
        let z = Zone::new("z", vec![s]).unwrap();
        let building = Building::new("b", vec![z]).unwrap();
        let surface_index = SurfaceIndex::new(&building);
        let ctx = SimContext::new(&building, &surface_index);
        let mut bus = Bus::new();

        let module_cfg: EnergyModuleConfig = Default::default();
        let mut module = EnergyModule::new(module_cfg);
        let err = module.step(&ctx, &mut bus).unwrap_err();
        assert!(err.to_string().contains("not initialized"));
    }

    #[test]
    fn test_energy_module_consumes_shortwave_transmitted_per_zone() {
        // Two zones with no envelope loss and no coupling; shortwave in z0 should only heat z0.
        let s0 = Solid::from_box(1.0, 1.0, 1.0, None, "s0").unwrap();
        let s1 = Solid::from_box(1.0, 1.0, 1.0, Some((3.0, 0.0, 0.0)), "s1").unwrap();
        let z0 = Zone::new("z0", vec![s0]).unwrap();
        let z1 = Zone::new("z1", vec![s1]).unwrap();
        let building = Building::new("b", vec![z0, z1]).unwrap();

        let surface_index = SurfaceIndex::new(&building);
        let ctx = SimContext::new(&building, &surface_index);
        let mut bus = Bus::new();

        let mut thermal = ThermalConfig::new();
        thermal.default_u_value = 0.0;
        thermal.infiltration_ach = 0.0;
        thermal.indoor_temperature = 20.0;
        thermal.outdoor_temperature = 20.0;
        thermal.thermal_capacity_j_per_m3_k = 50_000.0;

        let hvac = HvacIdealLoads::with_setpoints(-1e9, 1e9);

        let mut module = EnergyModule::new(EnergyModuleConfig {
            thermal,
            hvac,
            dt_s: 3600.0,
            steady_state: false,
            model_kind: EnergyModelKind::AirOnly,
            material_library: None,
            default_envelope_capacity_j_per_m2_k: 0.0,
        });
        module.init(&ctx, &mut bus).unwrap();

        bus.put(OutdoorAirTemperatureC(20.0));

        let mut sw = ShortwaveTransmittedWPerZone::default();
        sw.watts_by_zone_uid
            .insert(module.zone_uids[0].clone(), 1000.0);
        bus.put(sw);

        module.step(&ctx, &mut bus).unwrap();

        let out = bus.get::<MultiZoneStepResult>().unwrap();
        assert!(out.zone_temperatures_c[0] > 20.0);
        assert!((out.zone_temperatures_c[1] - 20.0).abs() < 1e-6);
    }

    #[test]
    fn test_gains_by_zone_split_routes_absorbed_shortwave_by_exteriority() {
        // Two adjacent solids in the same zone => shared face is not exterior.
        let s0 = Solid::from_box(1.0, 1.0, 1.0, None, "s0").unwrap();
        let s1 = Solid::from_box(1.0, 1.0, 1.0, Some((1.0, 0.0, 0.0)), "s1").unwrap();
        let z = Zone::new("z", vec![s0, s1]).unwrap();
        let building = Building::new("b", vec![z]).unwrap();

        let surface_index = SurfaceIndex::new(&building);
        let ctx = SimContext::new(&building, &surface_index);
        let mut bus = Bus::new();

        let mut thermal = ThermalConfig::new();
        thermal.default_u_value = 0.0;
        thermal.infiltration_ach = 0.0;
        thermal.indoor_temperature = 20.0;
        thermal.outdoor_temperature = 20.0;

        let hvac = HvacIdealLoads::with_setpoints(-1e9, 1e9);
        let mut module = EnergyModule::new(EnergyModuleConfig {
            thermal,
            hvac,
            dt_s: 3600.0,
            steady_state: false,
            model_kind: EnergyModelKind::AirOnly,
            material_library: None,
            default_envelope_capacity_j_per_m2_k: 0.0,
        });
        module.init(&ctx, &mut bus).unwrap();

        let interior = surface_index
            .polygon_uid_by_path("z/s0/wall_1/poly_1")
            .unwrap()
            .clone();
        let exterior = surface_index
            .polygon_uid_by_path("z/s0/floor/floor")
            .unwrap()
            .clone();
        let boundaries = module.boundaries.as_ref().unwrap();
        assert!(!boundaries.is_exterior(&interior));
        assert!(boundaries.is_exterior(&exterior));

        let mut internal = InternalGainsWPerZone::default();
        internal
            .watts_by_zone_uid
            .insert(module.zone_uids[0].clone(), 50.0);
        bus.put(internal);

        let mut transmitted = ShortwaveTransmittedWPerZone::default();
        transmitted
            .watts_by_zone_uid
            .insert(module.zone_uids[0].clone(), 25.0);
        bus.put(transmitted);

        let mut absorbed = ShortwaveAbsorbedWPerPolygon::default();
        absorbed.watts_by_polygon_uid.insert(interior, 100.0);
        absorbed.watts_by_polygon_uid.insert(exterior, 200.0);
        bus.put(absorbed);

        let (air_gains, env_gains, internal_mass_sources) = module.gains_by_zone_split(&bus);
        assert_eq!(air_gains.len(), 1);
        assert_eq!(env_gains.len(), 1);
        assert_eq!(internal_mass_sources.len(), 1);
        assert!((air_gains[0] - (50.0 + 25.0 + 100.0)).abs() < 1e-12);
        assert!((env_gains[0] - 200.0).abs() < 1e-12);
        assert!((internal_mass_sources[0] - 0.0).abs() < 1e-12);
    }

    #[test]
    fn test_gains_by_zone_split_internal_total_with_zero_total_volume_returns_zero() {
        // Cover the defensive branch when `total_volume_m3` is effectively zero.
        let s0 = Solid::from_box(1.0, 1.0, 1.0, None, "s0").unwrap();
        let s1 = Solid::from_box(1.0, 1.0, 1.0, Some((3.0, 0.0, 0.0)), "s1").unwrap();
        let z0 = Zone::new("z0", vec![s0]).unwrap();
        let z1 = Zone::new("z1", vec![s1]).unwrap();
        let building = Building::new("b", vec![z0, z1]).unwrap();

        let surface_index = SurfaceIndex::new(&building);
        let ctx = SimContext::new(&building, &surface_index);
        let mut bus = Bus::new();

        let mut module = EnergyModule::new(EnergyModuleConfig::default());
        module.init(&ctx, &mut bus).unwrap();
        assert!(module.total_volume_m3 > 1e-6);

        // Force the "else { 0.0 }" internal gains branch.
        module.total_volume_m3 = 0.0;
        bus.put(InternalGainsWTotal(123.0));

        let (air_gains, env_gains, internal_mass_sources) = module.gains_by_zone_split(&bus);
        assert_eq!(air_gains, vec![0.0, 0.0]);
        assert_eq!(env_gains, vec![0.0, 0.0]);
        assert_eq!(internal_mass_sources, vec![0.0, 0.0]);
    }

    #[test]
    fn test_energy_module_envelope_rc_steady_state_init_path() {
        let s = Solid::from_box(2.0, 2.0, 2.0, None, "s").unwrap();
        let z = Zone::new("z", vec![s]).unwrap();
        let building = Building::new("b", vec![z]).unwrap();

        let surface_index = SurfaceIndex::new(&building);
        let ctx = SimContext::new(&building, &surface_index);
        let mut bus = Bus::new();

        let cfg = EnergyModuleConfig {
            model_kind: EnergyModelKind::EnvelopeRc2R1C,
            steady_state: true,
            ..EnergyModuleConfig::default()
        };
        let mut module = EnergyModule::new(cfg);
        module.init(&ctx, &mut bus).unwrap();

        match module.model.as_ref().unwrap() {
            EnergyModel::EnvelopeRc(_) => {}
            _ => panic!("expected EnvelopeRc model"),
        }
    }

    #[test]
    fn test_energy_module_envelope_rc_has_thermal_lag() {
        // No gains, HVAC off, outdoor drops: air-only steady-state should jump to outdoor,
        // while RC envelope retains heat for at least one step.
        let s = Solid::from_box(3.0, 3.0, 3.0, None, "s").unwrap();
        let z = Zone::new("z", vec![s]).unwrap();
        let building = Building::new("b", vec![z]).unwrap();

        let surface_index = SurfaceIndex::new(&building);
        let ctx = SimContext::new(&building, &surface_index);

        let mut thermal = ThermalConfig::new();
        thermal.default_u_value = 1.0;
        thermal.infiltration_ach = 0.0;
        thermal.thermal_capacity_j_per_m3_k = 0.0; // isolate envelope capacity effect
        thermal.indoor_temperature = 20.0;
        thermal.outdoor_temperature = 20.0;

        let hvac = HvacIdealLoads::with_setpoints(-1e9, 1e9);

        let mut bus_air = Bus::new();
        let mut air_only = EnergyModule::new(EnergyModuleConfig {
            thermal: thermal.clone(),
            hvac: hvac.clone(),
            dt_s: 3600.0,
            steady_state: true,
            model_kind: EnergyModelKind::AirOnly,
            material_library: None,
            default_envelope_capacity_j_per_m2_k: 0.0,
        });
        air_only.init(&ctx, &mut bus_air).unwrap();
        bus_air.put(OutdoorAirTemperatureC(0.0));
        air_only.step(&ctx, &mut bus_air).unwrap();
        let out_air = bus_air.get::<MultiZoneStepResult>().unwrap();
        let t_air = out_air.zone_temperatures_c[0];

        let mut bus_rc = Bus::new();
        let mut rc = EnergyModule::new(EnergyModuleConfig {
            thermal,
            hvac,
            dt_s: 3600.0,
            steady_state: false,
            model_kind: EnergyModelKind::EnvelopeRc2R1C,
            material_library: None,
            default_envelope_capacity_j_per_m2_k: 200_000.0,
        });
        rc.init(&ctx, &mut bus_rc).unwrap();
        bus_rc.put(OutdoorAirTemperatureC(0.0));
        rc.step(&ctx, &mut bus_rc).unwrap();
        let out_rc = bus_rc.get::<MultiZoneStepResult>().unwrap();
        let t_rc = out_rc.zone_temperatures_c[0];

        assert!(t_air < 1e-6, "air-only steady-state should go to outdoor");
        assert!(t_rc > 0.1, "RC envelope should retain heat at first step");
    }

    #[test]
    fn test_energy_module_consumes_shortwave_absorbed_per_polygon() {
        // Two zones with no envelope loss and no coupling; absorbed shortwave on a polygon
        // in z0 should only heat z0 (currently treated as a zone gain).
        let s0 = Solid::from_box(1.0, 1.0, 1.0, None, "s0").unwrap();
        let s1 = Solid::from_box(1.0, 1.0, 1.0, Some((3.0, 0.0, 0.0)), "s1").unwrap();
        let z0 = Zone::new("z0", vec![s0]).unwrap();
        let z1 = Zone::new("z1", vec![s1]).unwrap();
        let building = Building::new("b", vec![z0, z1]).unwrap();

        let surface_index = SurfaceIndex::new(&building);
        let ctx = SimContext::new(&building, &surface_index);
        let mut bus = Bus::new();

        let mut thermal = ThermalConfig::new();
        thermal.default_u_value = 0.0;
        thermal.infiltration_ach = 0.0;
        thermal.indoor_temperature = 20.0;
        thermal.outdoor_temperature = 20.0;
        thermal.thermal_capacity_j_per_m3_k = 50_000.0;

        let hvac = HvacIdealLoads::with_setpoints(-1e9, 1e9);

        let mut module = EnergyModule::new(EnergyModuleConfig {
            thermal,
            hvac,
            dt_s: 3600.0,
            steady_state: false,
            model_kind: EnergyModelKind::AirOnly,
            material_library: None,
            default_envelope_capacity_j_per_m2_k: 0.0,
        });
        module.init(&ctx, &mut bus).unwrap();

        bus.put(OutdoorAirTemperatureC(20.0));

        let z0_uid = ctx.building.zones().first().unwrap().uid.clone();
        let polygon_uid_z0 = surface_index
            .surfaces
            .iter()
            .find(|s| s.zone_uid == z0_uid)
            .unwrap()
            .polygon_uid
            .clone();

        let mut abs = ShortwaveAbsorbedWPerPolygon::default();
        abs.watts_by_polygon_uid.insert(polygon_uid_z0, 1000.0);
        bus.put(abs);

        module.step(&ctx, &mut bus).unwrap();
        let out = bus.get::<MultiZoneStepResult>().unwrap();
        assert!(out.zone_temperatures_c[0] > 20.0);
        assert!((out.zone_temperatures_c[1] - 20.0).abs() < 1e-6);
    }
}
