use std::collections::{HashMap, HashSet};

use crate::Building;

use crate::sim::heat_transfer::{BoundaryCondition, FvmWallSolver, build_1d_mesh};
use crate::sim::index::SurfaceIndex;

use super::boundary::ThermalBoundaries;
use super::config::{InternalMassBoundary, ThermalConfig};
use super::hvac::{HvacIdealLoads, LumpedThermalModel};
use super::network::{MultiZoneAirModel, ThermalNetwork};
use super::schedule::InternalGainsProfile;
use super::solar_bridge::{
    GlazingTransmission, SolarGainConfig, SolarHourParams, TransmittedSolarSplit,
    compute_glazing_transmissions_with_materials, compute_solar_gains_per_zone_with_materials,
    compute_solar_gains_with_materials,
};
use super::view_factors::SurfaceHandle;
use super::weather::WeatherData;
use super::zone::calculate_heat_balance_with_boundaries;
use crate::UID;
use crate::sim::engine::FlatScene;
use crate::sim::lighting::solar::SolarPosition;

#[cfg(test)]
use super::zone::calculate_heat_balance;

/// Result of an annual energy simulation.
#[derive(Debug, Clone)]
pub struct AnnualResult {
    /// Hourly heating demand in W for each hour.
    pub hourly_heating: Vec<f64>,
    /// Hourly cooling demand in W for each hour.
    pub hourly_cooling: Vec<f64>,
    /// Total annual heating energy in kWh.
    pub annual_heating_kwh: f64,
    /// Total annual cooling energy in kWh.
    pub annual_cooling_kwh: f64,
    /// Peak heating demand in W.
    pub peak_heating: f64,
    /// Peak cooling demand in W.
    pub peak_cooling: f64,
    /// Monthly heating energy in kWh (12 values).
    pub monthly_heating_kwh: [f64; 12],
    /// Monthly cooling energy in kWh (12 values).
    pub monthly_cooling_kwh: [f64; 12],
    /// Minimum zone air temperature over the year [C].
    pub min_zone_temp_c: f64,
    /// Maximum zone air temperature over the year [C].
    pub max_zone_temp_c: f64,
    /// Hourly zone air temperatures [C] (building-level: zone 0 or conditioned zone).
    pub hourly_zone_temp_c: Vec<f64>,
    /// Number of zones in the simulation.
    pub num_zones: usize,
    /// Per-zone hourly temperatures [C], indexed as `[zone_idx][hour]`.
    pub per_zone_hourly_temp_c: Vec<Vec<f64>>,
    /// Per-zone minimum temperature [C] over the year.
    pub per_zone_min_temp_c: Vec<f64>,
    /// Per-zone maximum temperature [C] over the year.
    pub per_zone_max_temp_c: Vec<f64>,
}

/// Multi-zone transient simulation output (zone air node per `Zone`).
#[derive(Debug, Clone)]
pub struct MultiZoneAnnualResult {
    pub zone_uids: Vec<UID>,
    pub zone_names: Vec<String>,
    /// Temperatures per zone per hour, indexed as `[zone][hour]`.
    pub hourly_zone_temperatures_c: Vec<Vec<f64>>,
    /// Thermal HVAC heating power per zone per hour, indexed as `[zone][hour]`.
    pub hourly_zone_heating_w: Vec<Vec<f64>>,
    /// Thermal HVAC cooling power per zone per hour, indexed as `[zone][hour]`.
    pub hourly_zone_cooling_w: Vec<Vec<f64>>,
    /// Aggregated building-level totals (sums over zones).
    pub annual: AnnualResult,
}

/// Computes day of year from month and day.
fn day_of_year(month: u8, day: u8) -> u16 {
    const DAYS_BEFORE_MONTH: [u16; 12] = [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334];
    let m = (month as usize).saturating_sub(1).min(11);
    DAYS_BEFORE_MONTH[m] + day as u16
}

fn opaque_absorptance_for_path(
    path: &str,
    material_library: Option<&crate::sim::materials::MaterialLibrary>,
    default_absorptance: f64,
) -> f64 {
    let Some(lib) = material_library else {
        return default_absorptance;
    };
    let Some(mat) = lib.lookup(path) else {
        return default_absorptance;
    };
    let Some(opt) = mat.optical.as_ref() else {
        return default_absorptance;
    };

    let mut a = 0.0;
    for c in 0..3 {
        let absorb =
            (1.0 - opt.diffuse_reflectance[c] - opt.specular_reflectance[c] - opt.transmittance[c])
                .clamp(0.0, 1.0);
        a += absorb;
    }
    (a / 3.0).clamp(0.0, 1.0)
}

#[derive(Clone)]
struct FvmExteriorWall {
    zone_uid: UID,
    polygon_uid: UID,
    path: String,
    area_m2: f64,
    normal: crate::Vector,
    is_ground_coupled: bool,
    h_out_w_per_m2_k: f64,
    h_min_iso_interior: f64,
    solver: FvmWallSolver,
    /// Cached interior surface temperature from the last solver step.
    /// Uses proper half-cell interpolation (not cell centroid).
    cached_interior_surface_temp_c: f64,
}

#[derive(Clone)]
struct FvmInternalMassSurface {
    zone_uid: UID,
    /// Index of this mass surface in the config's `internal_mass_surfaces` list
    /// (within the zone). Used as `SurfaceHandle::InternalMass { index }`.
    mass_index: usize,
    area_m2: f64,
    boundary: InternalMassBoundary,
    cos_tilt: f64,
    h_min_iso: f64,
    solver: FvmWallSolver,
    /// Cached interior surface temperature from the last solver step.
    cached_interior_surface_temp_c: f64,
}

/// A non-FVM exterior surface (e.g. glazing) whose interior surface temperature
/// is estimated via steady-state heat balance for inclusion in MRT.
#[derive(Clone)]
struct SteadyStateExteriorSurface {
    #[allow(dead_code)]
    zone_uid: UID,
    polygon_uid: UID,
    area_m2: f64,
    u_value_w_per_m2_k: f64,
    h_in_w_per_m2_k: f64,
    is_glazing: bool,
    #[allow(dead_code)]
    cos_tilt: f64,
    /// Cached interior surface temperature [C], solved explicitly in the global system.
    cached_interior_surface_temp_c: f64,
}

#[allow(clippy::too_many_arguments)]
fn distribute_transmitted_solar_geometric_fvm(
    transmissions: &[GlazingTransmission],
    fvm_walls: &[FvmExteriorWall],
    wall_centroids_by_uid: &HashMap<UID, crate::Point>,
    scene: &FlatScene,
    scene_uid_by_index: &[UID],
    sun_dir: crate::Vector,
    use_beam_distribution: bool,
    interior_solar_absorptance: f64,
    route_diffuse_to_air: bool,
) -> (HashMap<UID, f64>, f64) {
    let mut wall_sources_w: HashMap<UID, f64> = HashMap::new();
    let mut air_residual_w = 0.0_f64;
    if transmissions.is_empty() || fvm_walls.is_empty() {
        return (wall_sources_w, air_residual_w);
    }

    let alpha = interior_solar_absorptance.clamp(0.0, 1.0);

    let mut area_by_uid: HashMap<UID, f64> = HashMap::new();
    let mut zone_by_uid: HashMap<UID, UID> = HashMap::new();
    let mut candidates_by_zone_uid: HashMap<UID, Vec<UID>> = HashMap::new();
    for w in fvm_walls {
        if w.area_m2 <= 0.0 {
            continue;
        }
        area_by_uid.insert(w.polygon_uid.clone(), w.area_m2);
        zone_by_uid.insert(w.polygon_uid.clone(), w.zone_uid.clone());
        candidates_by_zone_uid
            .entry(w.zone_uid.clone())
            .or_default()
            .push(w.polygon_uid.clone());
    }

    for src in transmissions {
        let mut diffuse_pool_w = src.diffuse_w.max(0.0);

        // Move source slightly into the zone to avoid immediate self-hit.
        let inward = crate::Vector::new(
            -src.outward_normal.dx,
            -src.outward_normal.dy,
            -src.outward_normal.dz,
        )
        .normalize()
        .unwrap_or(crate::Vector::new(0.0, 0.0, -1.0));
        let source_pt = src.centroid + inward * 1e-3;

        if use_beam_distribution && src.beam_w > 0.0 {
            let beam_w = src.beam_w.max(0.0);
            let ray_dir = crate::Vector::new(-sun_dir.dx, -sun_dir.dy, -sun_dir.dz);
            if let Ok(ray_dir) = ray_dir.normalize() {
                if let Some((hit_idx, _)) = scene.find_target_surface_global(source_pt, ray_dir) {
                    if let Some(hit_uid) = scene_uid_by_index.get(hit_idx)
                        && area_by_uid.contains_key(hit_uid)
                        && zone_by_uid.get(hit_uid) == Some(&src.zone_uid)
                    {
                        let absorbed = alpha * beam_w;
                        if absorbed > 0.0 {
                            *wall_sources_w.entry(hit_uid.clone()).or_insert(0.0) += absorbed;
                        }
                        diffuse_pool_w += (1.0 - alpha) * beam_w;
                    } else {
                        diffuse_pool_w += beam_w;
                    }
                } else {
                    diffuse_pool_w += beam_w;
                }
            } else {
                diffuse_pool_w += beam_w;
            }
        } else {
            diffuse_pool_w += src.beam_w.max(0.0);
        }

        if diffuse_pool_w <= 0.0 {
            continue;
        }
        if route_diffuse_to_air {
            air_residual_w += diffuse_pool_w;
            continue;
        }

        let Some(candidates) = candidates_by_zone_uid.get(&src.zone_uid) else {
            air_residual_w += diffuse_pool_w;
            continue;
        };

        let mut weights: Vec<(UID, f64)> = Vec::new();
        for uid in candidates {
            let Some(&target_pt) = wall_centroids_by_uid.get(uid) else {
                continue;
            };
            let to_target = crate::Vector::new(
                target_pt.x - source_pt.x,
                target_pt.y - source_pt.y,
                target_pt.z - source_pt.z,
            );
            let dist_m = to_target.length();
            if dist_m <= 1e-6 {
                continue;
            }
            let Ok(ray_dir) = to_target.normalize() else {
                continue;
            };
            let Some((hit_idx, _)) = scene.find_target_surface_global(source_pt, ray_dir) else {
                continue;
            };
            let Some(hit_uid) = scene_uid_by_index.get(hit_idx) else {
                continue;
            };
            if hit_uid != uid {
                continue;
            }

            let area_m2 = area_by_uid.get(uid).copied().unwrap_or(0.0);
            if area_m2 <= 0.0 {
                continue;
            }
            let weight = area_m2 / dist_m.max(0.5).powi(2);
            if weight > 0.0 {
                weights.push((uid.clone(), weight));
            }
        }

        if weights.is_empty() {
            // Fallback: zone-local area-proportional split.
            let area_sum: f64 = candidates
                .iter()
                .map(|uid| area_by_uid.get(uid).copied().unwrap_or(0.0))
                .sum();
            if area_sum > 0.0 {
                for uid in candidates {
                    let a = area_by_uid.get(uid).copied().unwrap_or(0.0);
                    if a <= 0.0 {
                        continue;
                    }
                    *wall_sources_w.entry(uid.clone()).or_insert(0.0) +=
                        diffuse_pool_w * (a / area_sum);
                }
            } else {
                air_residual_w += diffuse_pool_w;
            }
            continue;
        }

        let w_sum: f64 = weights.iter().map(|(_, w)| *w).sum();
        if w_sum > 0.0 {
            for (uid, w) in weights {
                *wall_sources_w.entry(uid).or_insert(0.0) += diffuse_pool_w * (w / w_sum);
            }
        } else {
            air_residual_w += diffuse_pool_w;
        }
    }

    (wall_sources_w, air_residual_w)
}

fn is_ground_coupled_exterior_surface(
    config: &ThermalConfig,
    path: &str,
    normal: &crate::Vector,
) -> bool {
    if config.ground_temperature_c.is_none() {
        return false;
    }
    if !config
        .ground_surface_patterns
        .iter()
        .any(|p| path.contains(p.as_str()))
    {
        return false;
    }
    // Only treat downward-facing surfaces as ground-coupled.
    normal.dz <= -0.5
}

fn collect_fvm_exterior_walls(
    building: &Building,
    config: &ThermalConfig,
    index: &SurfaceIndex,
    boundaries: &ThermalBoundaries,
    solar_config: Option<&SolarGainConfig>,
) -> (Vec<FvmExteriorWall>, HashSet<UID>) {
    if !config.use_fvm_walls {
        return (vec![], HashSet::new());
    }

    let mut walls = Vec::new();
    let mut skip: HashSet<UID> = HashSet::new();

    for s in &index.surfaces {
        if !boundaries.is_exterior(&s.polygon_uid) {
            continue;
        }
        if s.area_m2 <= 0.0 {
            continue;
        }
        if looks_like_glazing(&s.path, config, solar_config) {
            continue;
        }
        if config.has_u_value_override_for_surface(&s.polygon_uid, &s.path) {
            continue;
        }

        let Some(construction) = config.resolve_construction(&s.path) else {
            continue;
        };

        let Some(poly) = building.get_polygon(&s.path) else {
            continue;
        };

        let is_ground_coupled = is_ground_coupled_exterior_surface(config, &s.path, &poly.vn);

        // Interior heat transfer coefficient.
        //
        // The ISO 6946 inside surface resistance `R_si` represents a *combined*
        // convection+radiation surface resistance. Since we do not yet model an explicit
        // interior longwave enclosure balance, we treat `1/R_si` as a reasonable first-order
        // combined coefficient for the FVM wall boundary.
        //
        // If the user provides an interior coefficient, treat it as a *minimum* (to avoid
        // unintentionally making the envelope artificially insulating).
        let h_si_iso = if construction.r_si > 0.0 {
            1.0 / construction.r_si
        } else {
            1.0 / 0.13
        }
        .max(1e-9);
        let h_out = if construction.r_se > 0.0 {
            1.0 / construction.r_se
        } else {
            1.0 / 0.04
        }
        .max(1e-9);

        let mesh = build_1d_mesh(construction, s.area_m2);
        let solver = FvmWallSolver::new(mesh, config.indoor_temperature);

        let init_temp = solver.interior_surface_temp();
        walls.push(FvmExteriorWall {
            zone_uid: s.zone_uid.clone(),
            polygon_uid: s.polygon_uid.clone(),
            path: s.path.clone(),
            area_m2: s.area_m2,
            normal: poly.vn,
            is_ground_coupled,
            h_out_w_per_m2_k: h_out,
            h_min_iso_interior: h_si_iso,
            solver,
            cached_interior_surface_temp_c: init_temp,
        });
        skip.insert(s.polygon_uid.clone());
    }

    (walls, skip)
}

fn collect_internal_mass_surfaces(
    building: &Building,
    config: &ThermalConfig,
) -> Vec<FvmInternalMassSurface> {
    if config.internal_mass_surfaces.is_empty() {
        return vec![];
    }

    let mut out = Vec::new();
    let mut mass_index_counter = 0_usize;
    for zone in building.zones() {
        for m in &config.internal_mass_surfaces {
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
            let this_mass_index = mass_index_counter;
            mass_index_counter += 1;

            let mesh = build_1d_mesh(&m.construction, m.area_m2);
            let solver = FvmWallSolver::new(mesh, config.indoor_temperature);
            let init_temp_in = solver.interior_surface_temp();
            out.push(FvmInternalMassSurface {
                zone_uid: zone.uid.clone(),
                mass_index: this_mass_index,
                area_m2: m.area_m2,
                boundary: m.boundary,
                cos_tilt: m.cos_tilt,
                h_min_iso: h_si_iso,
                solver,
                cached_interior_surface_temp_c: init_temp_in,
            });
        }
    }
    out
}

fn collect_steady_state_exterior_surfaces(
    building: &Building,
    config: &ThermalConfig,
    index: &SurfaceIndex,
    boundaries: &ThermalBoundaries,
    fvm_skip_polygons: &HashSet<UID>,
    solar_config: Option<&SolarGainConfig>,
) -> Vec<SteadyStateExteriorSurface> {
    let mut out = Vec::new();
    for s in &index.surfaces {
        if !boundaries.is_exterior(&s.polygon_uid) {
            continue;
        }
        if s.area_m2 <= 0.0 {
            continue;
        }
        if fvm_skip_polygons.contains(&s.polygon_uid) {
            continue;
        }
        // Skip ground-coupled surfaces.
        if let Some(poly) = building.get_polygon(&s.path)
            && is_ground_coupled_exterior_surface(config, &s.path, &poly.vn)
        {
            continue;
        }
        let u = config.resolve_u_value_for_surface(&s.polygon_uid, &s.path);
        if u <= 0.0 {
            continue;
        }
        let h_si_iso = 1.0 / 0.13_f64;
        let h_in = if config.interior_heat_transfer_coeff_w_per_m2_k > 0.0 {
            config.interior_heat_transfer_coeff_w_per_m2_k.max(h_si_iso)
        } else {
            h_si_iso
        };
        let cos_tilt = building
            .get_polygon(&s.path)
            .map(|p| p.vn.dz)
            .unwrap_or(0.0);
        let is_glazing = looks_like_glazing(&s.path, config, solar_config);
        let denom = u + h_in;
        let init_temp = if denom > 0.0 {
            config.indoor_temperature
                - (u / denom) * (config.indoor_temperature - config.outdoor_temperature)
        } else {
            config.indoor_temperature
        };
        out.push(SteadyStateExteriorSurface {
            zone_uid: s.zone_uid.clone(),
            polygon_uid: s.polygon_uid.clone(),
            area_m2: s.area_m2,
            u_value_w_per_m2_k: u,
            h_in_w_per_m2_k: h_in,
            is_glazing,
            cos_tilt,
            cached_interior_surface_temp_c: init_temp,
        });
    }
    out
}

#[allow(clippy::too_many_arguments)]
fn step_fvm_exterior_walls_fill_gains_by_zone_uid<F: Fn(&UID) -> f64, G: Fn(&UID) -> f64>(
    walls: &mut [FvmExteriorWall],
    config: &ThermalConfig,
    solar_config: Option<&SolarGainConfig>,
    params: Option<&SolarHourParams>,
    interior_surface_sources_w_by_zone_uid: Option<&std::collections::HashMap<UID, f64>>,
    interior_surface_sources_w_by_polygon_uid: Option<&HashMap<UID, f64>>,
    floor_beam_sources_w_by_zone_uid: Option<&std::collections::HashMap<UID, f64>>,
    outdoor_temp_c: f64,
    wind_speed: f64,
    zone_air_temp_for: F,
    zone_radiant_temp_for: G,
    dt_s: f64,
    gains_out: &mut std::collections::HashMap<UID, f64>,
    per_surface_mrt: Option<&std::collections::HashMap<SurfaceHandle, f64>>,
    h_r_uniform: f64,
    _convective_only_air_gain: bool,
) {
    gains_out.clear();
    if walls.is_empty() {
        return;
    }

    let mut sun_above = false;
    let mut sun_dir = crate::Vector::new(0.0, 0.0, 1.0);
    if let Some(p) = params {
        let solar_pos = SolarPosition::calculate_from_local_time(
            p.latitude,
            p.longitude,
            p.timezone,
            p.day_of_year,
            p.local_time_hours,
        );
        sun_above = solar_pos.is_above_horizon();
        sun_dir = solar_pos.to_direction();
    }

    fn isotropic_sky_view_factor_from_nz(n_z: f64) -> f64 {
        (0.5 * (1.0 + n_z)).clamp(0.0, 1.0)
    }

    let mut target_area_by_zone_uid: std::collections::HashMap<UID, f64> =
        std::collections::HashMap::new();
    if let Some(src) = interior_surface_sources_w_by_zone_uid {
        let mut total_area_by_zone_uid: std::collections::HashMap<UID, f64> =
            std::collections::HashMap::new();
        for w in walls.iter() {
            if w.area_m2 <= 0.0 {
                continue;
            }
            *total_area_by_zone_uid
                .entry(w.zone_uid.clone())
                .or_insert(0.0) += w.area_m2;
        }
        for (zone_uid, w_total) in src {
            if *w_total == 0.0 {
                continue;
            }
            let total_area = total_area_by_zone_uid.get(zone_uid).copied().unwrap_or(0.0);
            if total_area > 0.0 {
                target_area_by_zone_uid.insert(zone_uid.clone(), total_area);
            }
        }
    }

    // Floor beam solar: distributed only to floor FVM walls (normal.dz <= -0.5).
    let mut floor_beam_area_by_zone_uid: std::collections::HashMap<UID, f64> =
        std::collections::HashMap::new();
    if let Some(src) = floor_beam_sources_w_by_zone_uid {
        for w in walls.iter() {
            if w.area_m2 <= 0.0 || w.normal.dz > -0.5 {
                continue;
            }
            if src.get(&w.zone_uid).copied().unwrap_or(0.0) != 0.0 {
                *floor_beam_area_by_zone_uid
                    .entry(w.zone_uid.clone())
                    .or_insert(0.0) += w.area_m2;
            }
        }
    }

    for w in walls {
        let t_air = zone_air_temp_for(&w.zone_uid);
        // Dynamic interior convection coefficient.
        let cos_tilt = w.normal.dz;
        let t_surf_prev = w.solver.interior_surface_temp();
        let h_in_base = super::convection::interior_convection_h(
            &config.interior_convection_model,
            t_surf_prev - t_air,
            cos_tilt,
            w.h_min_iso_interior,
        )
        .max(1e-9);

        const SIGMA: f64 = 5.670_374_419e-8; // Stefan–Boltzmann (W/m²/K⁴)

        // Separate convective (zone air coupling) from total (FVM BC) coefficient.
        let (h_in_conv, h_in_total, t_eff) = if let Some(mrts) = per_surface_mrt {
            // View-factor path: TARP convection + uniform linearized radiation.
            let h_conv = h_in_base;
            let eps = config.interior_emissivity;
            let h_rad = eps * h_r_uniform;
            let t_mrt = mrts
                .get(&SurfaceHandle::Polygon(w.polygon_uid.clone()))
                .copied()
                .unwrap_or(t_air);
            let h_total = h_conv + h_rad;
            let t_eff = if h_total > 0.0 {
                (h_conv * t_air + h_rad * t_mrt) / h_total
            } else {
                t_air
            };
            (h_conv, h_total, t_eff)
        } else if config.use_interior_radiative_exchange {
            match config.interior_convection_model {
                super::convection::InteriorConvectionModel::Tarp => {
                    // TARP is convection-only; use it as the full coupling
                    // to air. Radiative exchange requires a proper view-factor
                    // model to avoid energy loss through incomplete MRT feedback.
                    (h_in_base, h_in_base, t_air)
                }
                super::convection::InteriorConvectionModel::Fixed(_) => {
                    let f_rad = config.interior_radiation_fraction.clamp(0.0, 1.0);
                    let h_rad = h_in_base * f_rad;
                    let h_conv = (h_in_base - h_rad).max(0.0);
                    let t_rad = zone_radiant_temp_for(&w.zone_uid);
                    let t_eff = if h_in_base > 0.0 {
                        (h_conv * t_air + h_rad * t_rad) / h_in_base
                    } else {
                        t_air
                    };
                    (h_conv, h_in_base, t_eff)
                }
            }
        } else {
            (h_in_base, h_in_base, t_air)
        };

        // Dynamic exterior convection coefficient.
        let t_surf_ext_prev = w.solver.exterior_surface_temp();
        let h_out_base = if let (Some(sc), Some(p)) = (solar_config, params) {
            let mut h = if sc.use_wind_speed_for_h_out {
                let v = p.wind_speed.max(0.0);
                sc.h_out_base_w_per_m2_k + sc.h_out_wind_coeff_w_per_m2_k_per_m_s * v
            } else {
                w.h_out_w_per_m2_k
            };
            if sc.h_out_tilt_scale != 0.0 {
                h *= 1.0 + sc.h_out_tilt_scale * w.normal.dz.abs();
            }
            h.max(1e-9)
        } else {
            w.h_out_w_per_m2_k.max(1e-9)
        };
        let h_out = super::convection::exterior_convection_h(
            &config.exterior_convection_model,
            h_out_base,
            t_surf_ext_prev - outdoor_temp_c,
            cos_tilt,
            wind_speed,
        )
        .max(1e-9);

        let mut heat_flux_sw_w_per_m2 = 0.0;

        if !w.is_ground_coupled
            && let (Some(sc), Some(p)) = (solar_config, params)
            && sc.include_exterior_opaque_absorption
        {
            let normal = w.normal;
            let sky_view = isotropic_sky_view_factor_from_nz(normal.dz.clamp(-1.0, 1.0));
            let mut incident = p.diffuse_horizontal_irradiance.max(0.0) * sky_view;
            if sun_above && p.direct_normal_irradiance > 0.0 {
                let cos_incidence = sun_dir.dot(&normal).max(0.0);
                incident += p.direct_normal_irradiance.max(0.0) * cos_incidence;
            }

            if incident > 0.0 {
                let a = opaque_absorptance_for_path(
                    &w.path,
                    config.material_library.as_ref(),
                    sc.default_opaque_absorptance,
                );
                if a > 0.0 {
                    heat_flux_sw_w_per_m2 = incident * a;
                }
            }
        }

        // Exterior longwave radiation exchange with sky/ground
        let mut heat_flux_net_w_per_m2 = heat_flux_sw_w_per_m2;
        if !w.is_ground_coupled
            && let (Some(sc), Some(p)) = (solar_config, params)
            && sc.include_exterior_longwave_exchange
        {
            let eps = sc.exterior_opaque_emissivity.clamp(0.0, 1.0);
            if eps > 0.0 {
                let normal = w.normal;
                let n_z = normal.dz.clamp(-1.0, 1.0);
                let sky_view = 0.5 * (1.0 + n_z);
                let ground_view = 1.0 - sky_view;
                let t_air_k = (p.outdoor_air_temperature_c + 273.15).max(1.0);
                let l_sky = p.horizontal_infrared_radiation.max(0.0);
                let eps_g = sc.ground_emissivity.clamp(0.0, 1.0);
                let l_ground = eps_g * SIGMA * t_air_k.powi(4);
                let incoming = sky_view * l_sky + ground_view * l_ground;
                let outgoing = SIGMA * t_air_k.powi(4);
                // q_lw < 0 when sky is cold → net heat loss from surface
                heat_flux_net_w_per_m2 += eps * (incoming - outgoing);
            }
        }

        let mut interior_source_flux_w_per_m2 = 0.0;
        if let Some(src) = interior_surface_sources_w_by_polygon_uid
            && let Some(&w_src) = src.get(&w.polygon_uid)
            && w_src != 0.0
            && w.area_m2 > 0.0
        {
            interior_source_flux_w_per_m2 += w_src / w.area_m2;
        }
        if let Some(src) = interior_surface_sources_w_by_zone_uid
            && let Some(&w_total) = src.get(&w.zone_uid)
            && w_total != 0.0
            && let Some(&target_area) = target_area_by_zone_uid.get(&w.zone_uid)
            && target_area > 0.0
        {
            interior_source_flux_w_per_m2 = w_total / target_area;
        }
        // Floor beam solar: additional flux only for floor walls (normal.dz <= -0.5).
        if w.normal.dz <= -0.5
            && let Some(src) = floor_beam_sources_w_by_zone_uid
            && let Some(&beam_w) = src.get(&w.zone_uid)
            && beam_w != 0.0
            && let Some(&floor_area) = floor_beam_area_by_zone_uid.get(&w.zone_uid)
            && floor_area > 0.0
        {
            interior_source_flux_w_per_m2 += beam_w / floor_area;
        }

        let bc_interior = if interior_source_flux_w_per_m2 != 0.0 {
            BoundaryCondition::ConvectiveWithFluxToDomain {
                h: h_in_total,
                t_fluid: t_eff,
                heat_flux: interior_source_flux_w_per_m2,
            }
        } else {
            BoundaryCondition::Convective {
                h: h_in_total,
                t_fluid: t_eff,
            }
        };

        let bc_exterior = if w.is_ground_coupled {
            BoundaryCondition::Convective {
                h: h_out,
                t_fluid: config.ground_temperature_c.unwrap_or(outdoor_temp_c),
            }
        } else if heat_flux_net_w_per_m2 != 0.0 {
            BoundaryCondition::ConvectiveWithFlux {
                h: h_out,
                t_fluid: outdoor_temp_c,
                heat_flux: heat_flux_net_w_per_m2,
            }
        } else {
            BoundaryCondition::Convective {
                h: h_out,
                t_fluid: outdoor_temp_c,
            }
        };

        w.solver.step(dt_s, &bc_exterior, &bc_interior, &[]);

        // Cache proper surface temperature for MRT computation next substep.
        w.cached_interior_surface_temp_c = w.solver.interior_surface_temperature(&bc_interior);

        // Report heat transfer to zone air.
        let t_surf = w.cached_interior_surface_temp_c;
        // Air node receives convective exchange only.
        // In view-factor mode, radiative exchange is surface-to-surface.
        let q_w_per_m2 = h_in_conv * (t_surf - t_air);
        if q_w_per_m2 != 0.0 {
            *gains_out.entry(w.zone_uid.clone()).or_insert(0.0) += q_w_per_m2 * w.area_m2;
        }
    }
}

fn compute_exterior_opaque_sol_air_gain_total_w(
    building: &Building,
    base_config: &ThermalConfig,
    index: &SurfaceIndex,
    boundaries: &ThermalBoundaries,
    params: &SolarHourParams,
    solar_config: &SolarGainConfig,
    skip_polygon_uids: Option<&HashSet<UID>>,
) -> f64 {
    const SIGMA: f64 = 5.670_374_419e-8; // Stefan–Boltzmann (W/m²/K⁴)
    fn isotropic_sky_view_factor_from_nz(n_z: f64) -> f64 {
        (0.5 * (1.0 + n_z)).clamp(0.0, 1.0)
    }

    let h_out = if solar_config.use_wind_speed_for_h_out {
        let v = params.wind_speed.max(0.0);
        solar_config.h_out_base_w_per_m2_k + solar_config.h_out_wind_coeff_w_per_m2_k_per_m_s * v
    } else {
        solar_config.exterior_heat_transfer_coeff_w_per_m2_k
    };

    let solar_pos = SolarPosition::calculate_from_local_time(
        params.latitude,
        params.longitude,
        params.timezone,
        params.day_of_year,
        params.local_time_hours,
    );
    let sun_above = solar_pos.is_above_horizon();
    let sun_dir = solar_pos.to_direction();

    let h_out = h_out.max(1e-9);

    let mut total = 0.0;
    for surface in &index.surfaces {
        if let Some(skip) = skip_polygon_uids
            && skip.contains(&surface.polygon_uid)
        {
            continue;
        }
        if !boundaries.is_exterior(&surface.polygon_uid) {
            continue;
        }
        if solar_config
            .resolve_shgc(&surface.path, base_config.material_library.as_ref())
            .is_some()
        {
            continue; // glazing handled via transmitted-to-zone path
        }

        let area = surface.area_m2;
        if area <= 0.0 {
            continue;
        }
        let Some(poly) = building.get_polygon(&surface.path) else {
            continue;
        };

        let normal = poly.vn;
        let n_z = normal.dz.clamp(-1.0, 1.0);
        let sky_view = isotropic_sky_view_factor_from_nz(n_z);
        let sky_view_lw = 0.5 * (1.0 + n_z);
        let ground_view_lw = 1.0 - sky_view_lw;

        let mut incident = params.diffuse_horizontal_irradiance.max(0.0) * sky_view;
        if sun_above && params.direct_normal_irradiance > 0.0 {
            let cos_incidence = sun_dir.dot(&normal).max(0.0);
            incident += params.direct_normal_irradiance.max(0.0) * cos_incidence;
        }

        let mut q_sw_w_per_m2 = 0.0;
        if solar_config.include_exterior_opaque_absorption && incident > 0.0 {
            let absorptance = opaque_absorptance_for_path(
                &surface.path,
                base_config.material_library.as_ref(),
                solar_config.default_opaque_absorptance,
            );
            if absorptance > 0.0 {
                q_sw_w_per_m2 += incident * absorptance;
            }
        }

        let mut h_conv = h_out;
        if solar_config.h_out_tilt_scale != 0.0 {
            h_conv *= 1.0 + solar_config.h_out_tilt_scale * n_z.abs();
            h_conv = h_conv.max(1e-9);
        }
        let mut q_lw_w_per_m2 = 0.0;
        let mut h_lw_total = h_conv;

        if solar_config.include_exterior_longwave_exchange {
            let eps = solar_config.exterior_opaque_emissivity.clamp(0.0, 1.0);
            if eps > 0.0 {
                let t_air_k = (params.outdoor_air_temperature_c + 273.15).max(1.0);
                let l_sky = params.horizontal_infrared_radiation.max(0.0);
                let eps_g = solar_config.ground_emissivity.clamp(0.0, 1.0);
                let l_ground = eps_g * SIGMA * t_air_k.powi(4);
                let incoming = sky_view_lw.max(0.0) * l_sky + ground_view_lw.max(0.0) * l_ground;
                let outgoing = SIGMA * t_air_k.powi(4);
                q_lw_w_per_m2 += eps * (incoming - outgoing);

                // Linearized radiative exchange coefficient at outdoor air temperature.
                // This helps avoid over-coupling longwave terms when using a single
                // sol-air-style gain path.
                let h_rad = (4.0 * eps * SIGMA * t_air_k.powi(3)).max(0.0);
                h_lw_total += h_rad;
            }
        }

        if q_sw_w_per_m2 == 0.0 && q_lw_w_per_m2 == 0.0 {
            continue;
        }

        let u = base_config.resolve_u_value_for_surface(&surface.polygon_uid, &surface.path);
        if !(u.is_finite() && u > 0.0) {
            continue;
        }

        // Sol-air style coupling: exterior fluxes shift the exterior surface temperature,
        // but only a fraction conducts inward instantaneously. Keep the historical shortwave
        // coupling behavior, and attenuate the longwave term by including an effective
        // radiative coefficient in the denominator.
        total += (u / h_conv) * q_sw_w_per_m2 * area + (u / h_lw_total) * q_lw_w_per_m2 * area;
    }
    total
}

fn compute_exterior_opaque_sol_air_gains_by_zone_w(
    building: &Building,
    base_config: &ThermalConfig,
    index: &SurfaceIndex,
    boundaries: &ThermalBoundaries,
    params: &SolarHourParams,
    solar_config: &SolarGainConfig,
    skip_polygon_uids: Option<&HashSet<UID>>,
) -> std::collections::HashMap<UID, f64> {
    const SIGMA: f64 = 5.670_374_419e-8;
    fn isotropic_sky_view_factor_from_nz(n_z: f64) -> f64 {
        (0.5 * (1.0 + n_z)).clamp(0.0, 1.0)
    }

    let h_out = if solar_config.use_wind_speed_for_h_out {
        let v = params.wind_speed.max(0.0);
        solar_config.h_out_base_w_per_m2_k + solar_config.h_out_wind_coeff_w_per_m2_k_per_m_s * v
    } else {
        solar_config.exterior_heat_transfer_coeff_w_per_m2_k
    };

    let solar_pos = SolarPosition::calculate_from_local_time(
        params.latitude,
        params.longitude,
        params.timezone,
        params.day_of_year,
        params.local_time_hours,
    );
    let sun_above = solar_pos.is_above_horizon();
    let sun_dir = solar_pos.to_direction();

    let h_out = h_out.max(1e-9);

    let mut gains: std::collections::HashMap<UID, f64> = std::collections::HashMap::new();
    for surface in &index.surfaces {
        if let Some(skip) = skip_polygon_uids
            && skip.contains(&surface.polygon_uid)
        {
            continue;
        }
        if !boundaries.is_exterior(&surface.polygon_uid) {
            continue;
        }
        if solar_config
            .resolve_shgc(&surface.path, base_config.material_library.as_ref())
            .is_some()
        {
            continue;
        }

        let area = surface.area_m2;
        if area <= 0.0 {
            continue;
        }
        let Some(poly) = building.get_polygon(&surface.path) else {
            continue;
        };

        let normal = poly.vn;
        let n_z = normal.dz.clamp(-1.0, 1.0);
        let sky_view = isotropic_sky_view_factor_from_nz(n_z);
        let sky_view_lw = 0.5 * (1.0 + n_z);
        let ground_view_lw = 1.0 - sky_view_lw;

        let mut incident = params.diffuse_horizontal_irradiance.max(0.0) * sky_view;
        if sun_above && params.direct_normal_irradiance > 0.0 {
            let cos_incidence = sun_dir.dot(&normal).max(0.0);
            incident += params.direct_normal_irradiance.max(0.0) * cos_incidence;
        }

        let mut q_sw_w_per_m2 = 0.0;
        if solar_config.include_exterior_opaque_absorption && incident > 0.0 {
            let absorptance = opaque_absorptance_for_path(
                &surface.path,
                base_config.material_library.as_ref(),
                solar_config.default_opaque_absorptance,
            );
            if absorptance > 0.0 {
                q_sw_w_per_m2 += incident * absorptance;
            }
        }

        let mut h_conv = h_out;
        if solar_config.h_out_tilt_scale != 0.0 {
            h_conv *= 1.0 + solar_config.h_out_tilt_scale * n_z.abs();
            h_conv = h_conv.max(1e-9);
        }
        let mut q_lw_w_per_m2 = 0.0;
        let mut h_lw_total = h_conv;

        if solar_config.include_exterior_longwave_exchange {
            let eps = solar_config.exterior_opaque_emissivity.clamp(0.0, 1.0);
            if eps > 0.0 {
                let t_air_k = (params.outdoor_air_temperature_c + 273.15).max(1.0);
                let l_sky = params.horizontal_infrared_radiation.max(0.0);
                let eps_g = solar_config.ground_emissivity.clamp(0.0, 1.0);
                let l_ground = eps_g * SIGMA * t_air_k.powi(4);
                let incoming = sky_view_lw.max(0.0) * l_sky + ground_view_lw.max(0.0) * l_ground;
                let outgoing = SIGMA * t_air_k.powi(4);
                q_lw_w_per_m2 += eps * (incoming - outgoing);

                let h_rad = (4.0 * eps * SIGMA * t_air_k.powi(3)).max(0.0);
                h_lw_total += h_rad;
            }
        }

        if q_sw_w_per_m2 == 0.0 && q_lw_w_per_m2 == 0.0 {
            continue;
        }

        let u = base_config.resolve_u_value_for_surface(&surface.polygon_uid, &surface.path);
        if !(u.is_finite() && u > 0.0) {
            continue;
        }

        let q = (u / h_conv) * q_sw_w_per_m2 * area + (u / h_lw_total) * q_lw_w_per_m2 * area;
        if q != 0.0 {
            *gains.entry(surface.zone_uid.clone()).or_insert(0.0) += q;
        }
    }
    gains
}

fn looks_like_glazing(
    path: &str,
    base_config: &ThermalConfig,
    solar_config: Option<&SolarGainConfig>,
) -> bool {
    if let Some(sc) = solar_config
        && sc
            .resolve_shgc(path, base_config.material_library.as_ref())
            .is_some()
    {
        return true;
    }
    if let Some(lib) = base_config.material_library.as_ref()
        && let Some(mat) = lib.lookup(path)
        && mat.is_glazing
    {
        return true;
    }
    let p = path.to_ascii_lowercase();
    p.contains("window") || p.contains("glazing") || p.contains("glass")
}

/// Runs an annual hourly energy simulation.
///
/// For each hour:
///   1. Set outdoor temperature from weather data
///   2. Calculate internal gains from schedule
///   3. Compute solar gains from DNI/DHI + sun geometry + window SHGC
///   4. Run steady-state heat balance
///   5. Accumulate heating/cooling demand
///
/// If `solar_config` is `None`, solar gains are zero for all hours.
pub fn run_annual_simulation(
    building: &Building,
    base_config: &ThermalConfig,
    weather: &WeatherData,
    gains_profile: Option<&InternalGainsProfile>,
    solar_config: Option<&SolarGainConfig>,
) -> AnnualResult {
    let index = SurfaceIndex::new(building);
    let boundaries = ThermalBoundaries::classify(building, &index);

    let num_hours = weather.num_hours();
    let mut hourly_heating = Vec::with_capacity(num_hours);
    let mut hourly_cooling = Vec::with_capacity(num_hours);

    let mut annual_heating = 0.0;
    let mut annual_cooling = 0.0;
    let mut peak_heating = 0.0_f64;
    let mut peak_cooling = 0.0_f64;
    let mut monthly_heating = [0.0; 12];
    let mut monthly_cooling = [0.0; 12];

    for (hour_idx, record) in weather.records.iter().enumerate() {
        let mut config = base_config.clone();
        config.outdoor_temperature = record.dry_bulb_temperature;

        // Internal gains
        if let Some(profile) = gains_profile {
            config.internal_gains = profile.gains_at(hour_idx);
        }

        // Solar gains from weather data + sun geometry + glazing properties (+ optional opaque sol-air coupling)
        config.solar_gains = match solar_config {
            Some(sc) => {
                let params = SolarHourParams {
                    outdoor_air_temperature_c: record.dry_bulb_temperature,
                    global_horizontal_irradiance: record.global_horizontal_radiation,
                    direct_normal_irradiance: record.direct_normal_radiation,
                    diffuse_horizontal_irradiance: record.diffuse_horizontal_radiation,
                    horizontal_infrared_radiation: record.horizontal_infrared_radiation,
                    wind_speed: record.wind_speed,
                    day_of_year: day_of_year(record.month, record.day),
                    local_time_hours: record.hour as f64 - 0.5,
                    latitude: weather.latitude,
                    longitude: weather.longitude,
                    timezone: weather.timezone,
                };
                let mut solar = compute_solar_gains_with_materials(
                    building,
                    &params,
                    sc,
                    config.material_library.as_ref(),
                )
                .total();
                if sc.include_exterior_opaque_absorption || sc.include_exterior_longwave_exchange {
                    solar += compute_exterior_opaque_sol_air_gain_total_w(
                        building,
                        &config,
                        &index,
                        &boundaries,
                        &params,
                        sc,
                        None,
                    );
                }
                solar
            }
            None => 0.0,
        };

        let result = calculate_heat_balance_with_boundaries(building, &config, &boundaries);

        hourly_heating.push(result.heating_demand);
        hourly_cooling.push(result.cooling_demand);

        annual_heating += result.heating_demand;
        annual_cooling += result.cooling_demand;
        peak_heating = peak_heating.max(result.heating_demand);
        peak_cooling = peak_cooling.max(result.cooling_demand);

        let month_idx = (record.month as usize).saturating_sub(1).min(11);
        monthly_heating[month_idx] += result.heating_demand;
        monthly_cooling[month_idx] += result.cooling_demand;
    }

    // Convert W*h to kWh
    let to_kwh = 1.0 / 1000.0;

    AnnualResult {
        hourly_heating,
        hourly_cooling,
        annual_heating_kwh: annual_heating * to_kwh,
        annual_cooling_kwh: annual_cooling * to_kwh,
        peak_heating,
        peak_cooling,
        monthly_heating_kwh: monthly_heating.map(|v| v * to_kwh),
        monthly_cooling_kwh: monthly_cooling.map(|v| v * to_kwh),
        min_zone_temp_c: 0.0,
        max_zone_temp_c: 0.0,
        hourly_zone_temp_c: Vec::new(),
        num_zones: 1,
        per_zone_hourly_temp_c: Vec::new(),
        per_zone_min_temp_c: Vec::new(),
        per_zone_max_temp_c: Vec::new(),
    }
}

/// Runs a transient annual simulation with thermal mass and HVAC.
///
/// Unlike `run_annual_simulation` which is steady-state, this models
/// the zone temperature state between time steps using a 1R1C model.
///
/// If `solar_config` is `None`, solar gains are zero for all hours.
pub fn run_transient_simulation(
    building: &Building,
    base_config: &ThermalConfig,
    weather: &WeatherData,
    hvac: &HvacIdealLoads,
    gains_profile: Option<&InternalGainsProfile>,
    solar_config: Option<&SolarGainConfig>,
) -> AnnualResult {
    run_transient_simulation_with_options(
        building,
        base_config,
        weather,
        hvac,
        gains_profile,
        solar_config,
        &TransientSimulationOptions::default(),
    )
}

#[derive(Debug, Clone)]
pub struct TransientSimulationOptions {
    /// Number of warmup hours to run before the reported simulation year.
    ///
    /// Warmup uses the first `warmup_hours` records of the provided weather data,
    /// then restarts reporting from hour 0 with the warmed-up state (similar in
    /// spirit to EnergyPlus warmup days).
    pub warmup_hours: usize,
    /// Number of internal substeps per EPW hour (default: 1).
    ///
    /// This is a numerical accuracy knob for transient simulations with explicit
    /// wall/internal-mass states. For heavyweight cases (e.g. BESTEST 900), using
    /// 6 (10-minute) or 12 (5-minute) substeps often improves stability and
    /// reduces timestep-induced bias compared to a single 1-hour step.
    pub substeps_per_hour: usize,
    /// Optional per-hour heating setpoint override [C]. Length must be >= num_hours.
    /// When `None`, uses `hvac.heating_setpoint`.
    pub hourly_heating_setpoint: Option<Vec<f64>>,
    /// Optional per-hour cooling setpoint override [C]. Length must be >= num_hours.
    /// When `None`, uses `hvac.cooling_setpoint`.
    pub hourly_cooling_setpoint: Option<Vec<f64>>,
    /// Optional per-hour infiltration ACH override.
    /// When `None`, uses `base_config.infiltration_ach`.
    pub hourly_infiltration_ach: Option<Vec<f64>>,
    /// Optional per-zone HVAC setpoints (indexed by zone index in sorted zone order).
    /// When `Some`, overrides the single `hvac` parameter for multi-zone simulations.
    /// Free-floating zones use extreme setpoints (heating -999, cooling 999).
    pub per_zone_hvac: Option<Vec<HvacIdealLoads>>,
}

impl Default for TransientSimulationOptions {
    fn default() -> Self {
        Self {
            warmup_hours: 0,
            substeps_per_hour: 1,
            hourly_heating_setpoint: None,
            hourly_cooling_setpoint: None,
            hourly_infiltration_ach: None,
            per_zone_hvac: None,
        }
    }
}

/// Like [`run_transient_simulation`], but supports an optional warmup period.
pub fn run_transient_simulation_with_options(
    building: &Building,
    base_config: &ThermalConfig,
    weather: &WeatherData,
    hvac: &HvacIdealLoads,
    gains_profile: Option<&InternalGainsProfile>,
    solar_config: Option<&SolarGainConfig>,
    options: &TransientSimulationOptions,
) -> AnnualResult {
    // Runtime policy:
    // - Keep only global FVM solve modes (VF off/on).
    // - Keep FVM wall conduction enabled.
    let mut enforced_config = base_config.clone();
    enforced_config.use_fvm_walls = true;
    let base_config = &enforced_config;

    let index = SurfaceIndex::new(building);
    let boundaries = ThermalBoundaries::classify(building, &index);
    let (mut fvm_walls, fvm_skip_polygons) =
        collect_fvm_exterior_walls(building, base_config, &index, &boundaries, solar_config);
    let has_fvm_walls = !fvm_walls.is_empty();

    let mut internal_mass_surfaces = collect_internal_mass_surfaces(building, base_config);
    let has_internal_mass = !internal_mass_surfaces.is_empty();

    let mut ss_exterior_surfaces = collect_steady_state_exterior_surfaces(
        building,
        base_config,
        &index,
        &boundaries,
        &fvm_skip_polygons,
        solar_config,
    );

    // ── View-factor radiation ──
    // When enabled, compute per-zone view factors for all interior-facing surfaces.
    // This replaces the simplified area-weighted MRT with per-surface MRT.
    use super::view_factors::{
        InternalMassInfo, SurfaceHandle, ViewFactorData, compute_building_view_factors,
    };
    let view_factor_data: Option<ViewFactorData> = if base_config.use_view_factor_radiation {
        let mass_infos: Vec<InternalMassInfo> = internal_mass_surfaces
            .iter()
            .map(|m| InternalMassInfo {
                index: m.mass_index,
                zone_uid: m.zone_uid.clone(),
                cos_tilt: m.cos_tilt,
                area_m2: m.area_m2,
            })
            .collect();
        Some(compute_building_view_factors(
            building,
            &index,
            &boundaries,
            &mass_infos,
            base_config.view_factor_rays_per_surface,
        ))
    } else {
        None
    };

    // Precompute eligible interior-surface areas per zone for distributing surface sources
    // (transmitted solar + radiant internal gains). This includes:
    // - interior faces of FVM exterior walls (eligible opaque exterior polygons),
    // - explicit internal mass slabs.
    //
    // Distributing across all eligible area reduces peak surface fluxes and better matches
    // the "radiant-to-surfaces" behavior in EnergyPlus.
    let mut interior_source_area_walls_by_zone_uid: std::collections::HashMap<UID, f64> =
        std::collections::HashMap::new();
    if base_config.distribute_transmitted_solar_to_fvm_walls {
        for w in &fvm_walls {
            if w.area_m2 <= 0.0 {
                continue;
            }
            *interior_source_area_walls_by_zone_uid
                .entry(w.zone_uid.clone())
                .or_insert(0.0) += w.area_m2;
        }
    }
    let mut interior_source_area_mass_by_zone_uid: std::collections::HashMap<UID, f64> =
        std::collections::HashMap::new();
    for m in &internal_mass_surfaces {
        if m.area_m2 <= 0.0 {
            continue;
        }
        *interior_source_area_mass_by_zone_uid
            .entry(m.zone_uid.clone())
            .or_insert(0.0) += m.area_m2;
    }

    // Geometry cache for interior shortwave distribution.
    let mut fvm_wall_centroids_by_uid: HashMap<UID, crate::Point> = HashMap::new();
    for w in &fvm_walls {
        if let Some(poly) = building.get_polygon(&w.path) {
            fvm_wall_centroids_by_uid.insert(w.polygon_uid.clone(), poly.centroid());
        }
    }
    let interior_solar_scene = if has_fvm_walls && base_config.use_surface_aware_solar_distribution
    {
        Some(FlatScene::new(building, 2.0, true))
    } else {
        None
    };
    let interior_solar_scene_uid_by_index: Vec<UID> = interior_solar_scene
        .as_ref()
        .map(|scene| scene.polygons.iter().map(|p| p.uid.clone()).collect())
        .unwrap_or_default();

    let num_hours = weather.num_hours();
    let mut hourly_heating = Vec::with_capacity(num_hours);
    let mut hourly_cooling = Vec::with_capacity(num_hours);

    // ── Zone index mapping ────────────────────────────────────────────────
    // building.zones() returns sorted by name, giving stable zone ordering.
    let zones = building.zones();
    let num_zones = zones.len();
    let zone_uid_to_idx: HashMap<UID, usize> = zones
        .iter()
        .enumerate()
        .map(|(i, z)| (z.uid.clone(), i))
        .collect();

    // Per-zone volumes and derived quantities.
    let zone_volumes: Vec<f64> = zones.iter().map(|z| z.volume()).collect();
    let volume: f64 = zone_volumes.iter().sum();

    let base_infiltration_ach = base_config.infiltration_ach;
    let zone_air_capacities: Vec<f64> = zone_volumes.iter().map(|v| 1.2 * 1005.0 * v).collect();
    let zone_infiltration_conds: Vec<f64> = zone_volumes
        .iter()
        .map(|v| 1.2 * 1005.0 * v * base_infiltration_ach / 3600.0)
        .collect();

    // Compute building-level UA breakdown from exterior surfaces.
    let mut ua_total = 0.0;
    let mut ua_ground = 0.0;
    // Per-zone glazing UA for the global solver.
    let mut zone_glazing_ua = vec![0.0_f64; num_zones];
    for s in &index.surfaces {
        if !boundaries.is_exterior(&s.polygon_uid) {
            continue;
        }
        if fvm_skip_polygons.contains(&s.polygon_uid) {
            continue;
        }
        let u = base_config.resolve_u_value_for_surface(&s.polygon_uid, &s.path);
        let ua = (u * s.area_m2).max(0.0);
        ua_total += ua;
        if looks_like_glazing(&s.path, base_config, solar_config) {
            let zi = zone_uid_to_idx.get(&s.zone_uid).copied().unwrap_or(0);
            zone_glazing_ua[zi] += ua;
        }

        if base_config.ground_temperature_c.is_some()
            && base_config
                .ground_surface_patterns
                .iter()
                .any(|p| s.path.contains(p.as_str()))
        {
            let Some(poly) = building.get_polygon(&s.path) else {
                continue;
            };
            if poly.vn.dz <= -0.5 {
                ua_ground += ua;
            }
        }
    }

    // For backward compat, single-zone aggregates.
    let infiltration_cond: f64 = zone_infiltration_conds.iter().sum();
    let air_capacity_j_per_k: f64 = zone_air_capacities.iter().sum();
    let thermal_capacity = air_capacity_j_per_k;

    let substeps_per_hour = options.substeps_per_hour.max(1);
    let dt_s = 3600.0 / substeps_per_hour as f64;
    let dt_h = dt_s / 3600.0;

    let mut model_1r1c = LumpedThermalModel::new(
        base_config.indoor_temperature,
        ua_total,
        infiltration_cond,
        thermal_capacity,
    );

    let mut annual_heating = 0.0;
    let mut annual_cooling = 0.0;
    let mut peak_heating = 0.0_f64;
    let mut peak_cooling = 0.0_f64;
    let mut monthly_heating = [0.0; 12];
    let mut monthly_cooling = [0.0; 12];
    let mut min_zone_temp = f64::MAX;
    let mut max_zone_temp = f64::MIN;
    let mut hourly_zone_temp: Vec<f64> = Vec::with_capacity(num_hours);
    let mut per_zone_hourly_temp: Vec<Vec<f64>> = vec![Vec::with_capacity(num_hours); num_zones];
    let mut per_zone_min_temp = vec![f64::MAX; num_zones];
    let mut per_zone_max_temp = vec![f64::MIN; num_zones];

    let warmup_hours = options.warmup_hours.min(num_hours);

    // ── Collect inter-zone partition walls ────────────────────────────────
    use super::boundary::ThermalSurfaceKind;
    let mut interzone_fvm_walls: Vec<FvmExteriorWall> = Vec::new();
    // Track which zone UID is on the "other" side of each inter-zone wall.
    let mut interzone_exterior_zone_uids: Vec<UID> = Vec::new();
    {
        let mut seen_pairs: HashSet<(UID, UID)> = HashSet::new();
        for (uid1, uid2) in &boundaries.facing_pairs {
            if boundaries.kind(uid1) != ThermalSurfaceKind::InterZoneInterface {
                continue;
            }
            // Avoid double-counting: pick one polygon per pair.
            // Use both orderings in the set since UID doesn't impl Ord.
            let key = (uid1.clone(), uid2.clone());
            let key_rev = (uid2.clone(), uid1.clone());
            if seen_pairs.contains(&key) || seen_pairs.contains(&key_rev) {
                continue;
            }
            seen_pairs.insert(key);
            // Use uid1 as the reference polygon.
            let Some(sref1) = index.surfaces.iter().find(|s| &s.polygon_uid == uid1) else {
                continue;
            };
            let Some(sref2) = index.surfaces.iter().find(|s| &s.polygon_uid == uid2) else {
                continue;
            };
            if sref1.area_m2 <= 0.0 {
                continue;
            }
            let Some(construction) = base_config.resolve_construction(&sref1.path) else {
                // Try uid2's path as fallback.
                let Some(c2) = base_config.resolve_construction(&sref2.path) else {
                    continue;
                };
                let mesh = build_1d_mesh(c2, sref1.area_m2);
                let solver = FvmWallSolver::new(mesh, base_config.indoor_temperature);
                let init_temp = solver.interior_surface_temp();
                let poly = building.get_polygon(&sref1.path);
                let normal = poly
                    .map(|p| p.vn)
                    .unwrap_or(crate::Vector::new(0.0, 0.0, 1.0));
                interzone_fvm_walls.push(FvmExteriorWall {
                    zone_uid: sref1.zone_uid.clone(),
                    polygon_uid: sref1.polygon_uid.clone(),
                    path: sref1.path.clone(),
                    area_m2: sref1.area_m2,
                    normal,
                    is_ground_coupled: false,
                    h_out_w_per_m2_k: 0.0,
                    h_min_iso_interior: 1.0 / 0.13,
                    solver,
                    cached_interior_surface_temp_c: init_temp,
                });
                interzone_exterior_zone_uids.push(sref2.zone_uid.clone());
                continue;
            };
            let mesh = build_1d_mesh(construction, sref1.area_m2);
            let solver = FvmWallSolver::new(mesh, base_config.indoor_temperature);
            let init_temp = solver.interior_surface_temp();
            let poly = building.get_polygon(&sref1.path);
            let normal = poly
                .map(|p| p.vn)
                .unwrap_or(crate::Vector::new(0.0, 0.0, 1.0));
            interzone_fvm_walls.push(FvmExteriorWall {
                zone_uid: sref1.zone_uid.clone(),
                polygon_uid: sref1.polygon_uid.clone(),
                path: sref1.path.clone(),
                area_m2: sref1.area_m2,
                normal,
                is_ground_coupled: false,
                h_out_w_per_m2_k: 0.0,
                h_min_iso_interior: if construction.r_si > 0.0 {
                    1.0 / construction.r_si
                } else {
                    1.0 / 0.13
                },
                solver,
                cached_interior_surface_temp_c: init_temp,
            });
            interzone_exterior_zone_uids.push(sref2.zone_uid.clone());
        }
    }
    let _has_interzone_walls = !interzone_fvm_walls.is_empty();

    // ── Global FVM solver initialization ────────────────────────────────
    use super::global_solve::{
        self, AirStepConditions, FvmWallInfo, GlobalTopology, RadiationConditions,
        SteadySurfaceInfo, WallStepConditions,
    };
    let mut global_ss_surface_map: Vec<usize> = Vec::new();
    let (mut global_topology, mut global_temps): (GlobalTopology, Vec<f64>) = {
        // Build wall info for topology
        let wall_infos: Vec<FvmWallInfo> = {
            let mut infos = Vec::new();
            for w in &fvm_walls {
                let zi = zone_uid_to_idx.get(&w.zone_uid).copied().unwrap_or(0);
                infos.push(FvmWallInfo {
                    solver: &w.solver,
                    zone_idx: zi,
                    area_m2: w.area_m2,
                    has_exterior_surface: false,
                    exterior_adiabatic: false,
                    exterior_zone_idx: None,
                });
            }
            for m in &internal_mass_surfaces {
                let zi = zone_uid_to_idx.get(&m.zone_uid).copied().unwrap_or(0);
                infos.push(FvmWallInfo {
                    solver: &m.solver,
                    zone_idx: zi,
                    area_m2: m.area_m2,
                    has_exterior_surface: matches!(m.boundary, InternalMassBoundary::TwoSided),
                    exterior_adiabatic: matches!(
                        m.boundary,
                        InternalMassBoundary::OneSidedAdiabatic
                    ),
                    exterior_zone_idx: None,
                });
            }
            // Inter-zone partition walls: TwoSided with exterior_zone_idx set.
            for (iw_idx, iw) in interzone_fvm_walls.iter().enumerate() {
                let zi = zone_uid_to_idx.get(&iw.zone_uid).copied().unwrap_or(0);
                let ext_zi = zone_uid_to_idx
                    .get(&interzone_exterior_zone_uids[iw_idx])
                    .copied()
                    .unwrap_or(0);
                infos.push(FvmWallInfo {
                    solver: &iw.solver,
                    zone_idx: zi,
                    area_m2: iw.area_m2,
                    has_exterior_surface: true,
                    exterior_adiabatic: false,
                    exterior_zone_idx: Some(ext_zi),
                });
            }
            infos
        };

        let mut steady_infos: Vec<SteadySurfaceInfo> = Vec::new();
        let mut explicit_glazing_ua_per_zone = vec![0.0_f64; num_zones];
        for (ss_idx, ss) in ss_exterior_surfaces.iter().enumerate() {
            if ss.area_m2 <= 0.0 || ss.u_value_w_per_m2_k <= 0.0 || ss.h_in_w_per_m2_k <= 0.0 {
                continue;
            }
            let zi = zone_uid_to_idx.get(&ss.zone_uid).copied().unwrap_or(0);
            if ss.is_glazing {
                explicit_glazing_ua_per_zone[zi] += ss.u_value_w_per_m2_k * ss.area_m2;
                zone_glazing_ua[zi] += ss.u_value_w_per_m2_k * ss.area_m2;
            }
            steady_infos.push(SteadySurfaceInfo {
                zone_idx: zi,
                area_m2: ss.area_m2,
                h_in_w_per_m2_k: ss.h_in_w_per_m2_k,
                u_to_out_w_per_m2_k: ss.u_value_w_per_m2_k,
            });
            global_ss_surface_map.push(ss_idx);
        }

        // Per-zone glazing UA residual (after subtracting explicitly represented glazing).
        let zone_glazing_ua_residual: Vec<f64> = zone_glazing_ua
            .iter()
            .zip(explicit_glazing_ua_per_zone.iter())
            .map(|(total, explicit)| (total - explicit).max(0.0))
            .collect();

        let topo = global_solve::build_topology_with_steady(
            &wall_infos,
            &steady_infos,
            num_zones,
            &zone_air_capacities,
            &zone_infiltration_conds,
            &zone_glazing_ua_residual,
        );

        // Extract initial temperatures
        let solvers: Vec<&FvmWallSolver> = {
            let mut s: Vec<&FvmWallSolver> = Vec::new();
            for w in &fvm_walls {
                s.push(&w.solver);
            }
            for m in &internal_mass_surfaces {
                s.push(&m.solver);
            }
            for iw in &interzone_fvm_walls {
                s.push(&iw.solver);
            }
            s
        };
        let air_temps_init: Vec<f64> = vec![model_1r1c.zone_temperature; num_zones];
        let mut steady_temps = Vec::with_capacity(global_ss_surface_map.len());
        for &ss_idx in &global_ss_surface_map {
            steady_temps.push(ss_exterior_surfaces[ss_idx].cached_interior_surface_temp_c);
        }
        let temps = global_solve::extract_temperatures_with_steady(
            &topo,
            &solvers,
            &air_temps_init,
            &steady_temps,
        );
        (topo, temps)
    };

    let mut simulate_hour = |hour_idx: usize,
                             record: &super::weather::HourlyRecord,
                             report: bool,
                             hour_hvac: &HvacIdealLoads,
                             hour_infiltration_ach: f64| {
        // Apply per-hour infiltration ACH override (per-zone).
        let hour_infiltration_cond = 1.2 * 1005.0 * volume * hour_infiltration_ach / 3600.0;
        model_1r1c.infiltration_conductance = hour_infiltration_cond;
        for air in &mut global_topology.air_nodes {
            let zv = zone_volumes.get(air.zone_idx).copied().unwrap_or(0.0);
            air.infiltration_k = 1.2 * 1005.0 * zv * hour_infiltration_ach / 3600.0;
        }
        let gains = gains_profile
            .map(|p| p.gains_at(hour_idx))
            .unwrap_or(base_config.internal_gains);
        let mut solar_params: Option<SolarHourParams> = None;
        let mut transmitted_split = TransmittedSolarSplit::default();
        let mut glazing_transmissions: Vec<GlazingTransmission> = Vec::new();
        let (solar_transmitted, solar_opaque_sol_air) = match solar_config {
            Some(sc) => {
                let params = SolarHourParams {
                    outdoor_air_temperature_c: record.dry_bulb_temperature,
                    global_horizontal_irradiance: record.global_horizontal_radiation,
                    direct_normal_irradiance: record.direct_normal_radiation,
                    diffuse_horizontal_irradiance: record.diffuse_horizontal_radiation,
                    horizontal_infrared_radiation: record.horizontal_infrared_radiation,
                    wind_speed: record.wind_speed,
                    day_of_year: day_of_year(record.month, record.day),
                    local_time_hours: record.hour as f64 - 0.5,
                    latitude: weather.latitude,
                    longitude: weather.longitude,
                    timezone: weather.timezone,
                };
                solar_params = Some(params);
                glazing_transmissions = compute_glazing_transmissions_with_materials(
                    building,
                    &params,
                    sc,
                    base_config.material_library.as_ref(),
                );
                let mut beam_w = 0.0_f64;
                let mut diffuse_w = 0.0_f64;
                for t in &glazing_transmissions {
                    beam_w += t.beam_w;
                    diffuse_w += t.diffuse_w;
                }
                transmitted_split = TransmittedSolarSplit { beam_w, diffuse_w };
                let transmitted = transmitted_split.total();
                let opaque = if sc.include_exterior_opaque_absorption
                    || sc.include_exterior_longwave_exchange
                {
                    compute_exterior_opaque_sol_air_gain_total_w(
                        building,
                        base_config,
                        &index,
                        &boundaries,
                        &params,
                        sc,
                        Some(&fvm_skip_polygons),
                    )
                } else {
                    0.0
                };
                (transmitted, opaque)
            }
            None => (0.0, 0.0),
        };
        let use_surface_sources = base_config.use_surface_aware_solar_distribution
            && (has_internal_mass || has_fvm_walls);

        let (gains_air_w, gains_surface_w) = if use_surface_sources {
            let f_surface = base_config.internal_gains_to_mass_fraction.clamp(0.0, 1.0);
            (gains * (1.0 - f_surface), gains * f_surface)
        } else {
            (gains, 0.0)
        };

        let (solar_transmitted_air_w, solar_transmitted_surface_w) = if use_surface_sources {
            let f_air = base_config
                .transmitted_solar_to_air_fraction
                .clamp(0.0, 1.0);
            (solar_transmitted * f_air, solar_transmitted * (1.0 - f_air))
        } else {
            (solar_transmitted, 0.0)
        };

        let mut solar_total_air_w = solar_transmitted_air_w + solar_opaque_sol_air;

        let mut interior_sources_walls_w_by_zone_uid: Option<std::collections::HashMap<UID, f64>> =
            None;
        let mut interior_sources_walls_w_by_polygon_uid: Option<HashMap<UID, f64>> = None;
        let mut interior_sources_mass_w_by_zone_uid: Option<std::collections::HashMap<UID, f64>> =
            None;
        let mut floor_beam_w_by_zone_uid: Option<std::collections::HashMap<UID, f64>> = None;
        if use_surface_sources {
            let f_air = base_config
                .transmitted_solar_to_air_fraction
                .clamp(0.0, 1.0);

            // Determine the zone UID for surface source distribution.
            let zone_uid_opt = internal_mass_surfaces
                .first()
                .map(|m| m.zone_uid.clone())
                .or_else(|| fvm_walls.first().map(|w| w.zone_uid.clone()));

            if has_internal_mass {
                // Existing path: internal mass surfaces present.
                let (solar_beam_to_mass_w, solar_diffuse_surface_w) = if base_config
                    .use_beam_solar_distribution
                    && base_config.distribute_transmitted_solar_to_fvm_walls
                {
                    (
                        transmitted_split.beam_w * (1.0 - f_air),
                        transmitted_split.diffuse_w * (1.0 - f_air),
                    )
                } else {
                    (0.0, solar_transmitted_surface_w)
                };

                let w_area_split = gains_surface_w + solar_diffuse_surface_w;
                let w_total = w_area_split + solar_beam_to_mass_w;

                if w_total != 0.0
                    && let Some(z) = zone_uid_opt.clone()
                {
                    let mut walls_src = std::collections::HashMap::new();
                    let mut mass_src = std::collections::HashMap::new();
                    let a_walls = interior_source_area_walls_by_zone_uid
                        .get(&z)
                        .copied()
                        .unwrap_or(0.0);
                    let a_mass = interior_source_area_mass_by_zone_uid
                        .get(&z)
                        .copied()
                        .unwrap_or(0.0);
                    let a_total = a_walls + a_mass;
                    if a_total > 0.0 {
                        let wall_w = if a_walls > 0.0 {
                            w_area_split * (a_walls / a_total)
                        } else {
                            0.0
                        };
                        let mass_w = w_area_split * (a_mass / a_total) + solar_beam_to_mass_w;

                        if wall_w > 0.0 {
                            if base_config.fvm_wall_solar_to_air {
                                solar_total_air_w += wall_w;
                            } else {
                                walls_src.insert(z.clone(), wall_w);
                            }
                        }
                        if mass_w > 0.0 {
                            mass_src.insert(z.clone(), mass_w);
                        }
                    }
                    if !walls_src.is_empty() {
                        interior_sources_walls_w_by_zone_uid = Some(walls_src);
                    }
                    if !mass_src.is_empty() {
                        interior_sources_mass_w_by_zone_uid = Some(mass_src);
                    }
                }
            } else if has_fvm_walls {
                // No internal mass: distribute radiant internal gains area-proportionally
                // to FVM walls, and distribute transmitted solar geometrically using
                // glazing source positions + ray visibility.
                let mut walls_src_by_polygon_uid: HashMap<UID, f64> = HashMap::new();

                if gains_surface_w != 0.0 {
                    let area_sum: f64 = fvm_walls
                        .iter()
                        .filter(|w| w.area_m2 > 0.0)
                        .map(|w| w.area_m2)
                        .sum();
                    if area_sum > 0.0 {
                        for w in &fvm_walls {
                            if w.area_m2 <= 0.0 {
                                continue;
                            }
                            *walls_src_by_polygon_uid
                                .entry(w.polygon_uid.clone())
                                .or_insert(0.0) += gains_surface_w * (w.area_m2 / area_sum);
                        }
                    }
                }

                if solar_transmitted_surface_w > 0.0 {
                    let solar_pos = solar_params.as_ref().map(|p| {
                        SolarPosition::calculate_from_local_time(
                            p.latitude,
                            p.longitude,
                            p.timezone,
                            p.day_of_year,
                            p.local_time_hours,
                        )
                    });
                    if let (Some(sp), Some(scene)) = (solar_pos, interior_solar_scene.as_ref()) {
                        let (solar_by_polygon_uid, air_residual_w) =
                            distribute_transmitted_solar_geometric_fvm(
                                &glazing_transmissions,
                                &fvm_walls,
                                &fvm_wall_centroids_by_uid,
                                scene,
                                &interior_solar_scene_uid_by_index,
                                sp.to_direction(),
                                base_config.use_beam_solar_distribution,
                                base_config.interior_solar_absorptance,
                                base_config.fvm_wall_solar_to_air,
                            );
                        for (uid, w_src) in solar_by_polygon_uid {
                            *walls_src_by_polygon_uid.entry(uid).or_insert(0.0) += w_src;
                        }
                        if air_residual_w > 0.0 {
                            solar_total_air_w += air_residual_w;
                        }
                    } else if let Some(z) = zone_uid_opt.clone() {
                        // Fallback: preserve previous area/floor heuristic if geometry
                        // inputs are unavailable.
                        let alpha = base_config.interior_solar_absorptance;
                        let (solar_beam_to_floor_w, solar_diffuse_surface_w) =
                            if base_config.use_beam_solar_distribution {
                                let raw_beam = transmitted_split.beam_w * (1.0 - f_air);
                                let raw_diffuse = transmitted_split.diffuse_w * (1.0 - f_air);
                                let floor_absorbed = alpha * raw_beam;
                                let reflected_pool = (1.0 - alpha) * raw_beam;
                                (floor_absorbed, raw_diffuse + reflected_pool)
                            } else {
                                (0.0, solar_transmitted_surface_w)
                            };
                        if solar_diffuse_surface_w > 0.0 {
                            if base_config.fvm_wall_solar_to_air {
                                solar_total_air_w += solar_diffuse_surface_w;
                            } else {
                                let mut walls_src = std::collections::HashMap::new();
                                walls_src.insert(z.clone(), solar_diffuse_surface_w);
                                interior_sources_walls_w_by_zone_uid = Some(walls_src);
                            }
                        }
                        if solar_beam_to_floor_w > 0.0 {
                            let mut floor_src = std::collections::HashMap::new();
                            floor_src.insert(z, solar_beam_to_floor_w);
                            floor_beam_w_by_zone_uid = Some(floor_src);
                        }
                    }
                }

                if !walls_src_by_polygon_uid.is_empty() {
                    interior_sources_walls_w_by_polygon_uid = Some(walls_src_by_polygon_uid);
                }
            }
        }
        let q_ground = if let Some(tg) = base_config.ground_temperature_c {
            ua_ground * (tg - record.dry_bulb_temperature)
        } else {
            0.0
        };

        let mut hour_heating_wh = 0.0_f64;
        let mut hour_cooling_wh = 0.0_f64;
        let mut hour_peak_heating_w = 0.0_f64;
        let mut hour_peak_cooling_w = 0.0_f64;

        for _substep in 0..substeps_per_hour {
            let (heating_power, cooling_power) = {
                // ── Global simultaneous FVM solve ───────────────────────
                // All wall cells + surface nodes + air node assembled into one matrix.
                let topo = &global_topology;
                let model = &mut model_1r1c;

                // Build per-wall conditions
                let mut g_wall_conds: Vec<WallStepConditions> = Vec::new();

                // Helper: compute solar position for this hour
                let mut sun_above = false;
                let mut sun_dir = crate::Vector::new(0.0, 0.0, 1.0);
                if let Some(p) = solar_params.as_ref() {
                    let solar_pos = SolarPosition::calculate_from_local_time(
                        p.latitude,
                        p.longitude,
                        p.timezone,
                        p.day_of_year,
                        p.local_time_hours,
                    );
                    sun_above = solar_pos.is_above_horizon();
                    sun_dir = solar_pos.to_direction();
                }

                fn isotropic_sky_view_g(n_z: f64) -> f64 {
                    (0.5 * (1.0 + n_z)).clamp(0.0, 1.0)
                }

                const SIGMA_G: f64 = 5.670_374_419e-8;

                let t_air_prev = model.zone_temperature;
                let t_out = record.dry_bulb_temperature;

                // ── Compute interior surface sources per wall ──────────
                // Collect total interior surface area for solar distribution
                let mut floor_area = 0.0f64;
                for w in fvm_walls.iter() {
                    if w.normal.dz <= -0.5 {
                        floor_area += w.area_m2;
                    }
                }
                for m in internal_mass_surfaces.iter() {
                    if m.cos_tilt <= -0.5 {
                        floor_area += m.area_m2;
                    }
                }

                // FVM exterior walls
                for w in fvm_walls.iter() {
                    let cos_tilt = w.normal.dz;
                    let t_surf_prev = w.solver.interior_surface_temp();
                    let h_in_base = super::convection::interior_convection_h(
                        &base_config.interior_convection_model,
                        t_surf_prev - t_air_prev,
                        cos_tilt,
                        w.h_min_iso_interior,
                    )
                    .max(1e-9);

                    // Exterior convection
                    let t_surf_ext_prev = w.solver.exterior_surface_temp();
                    let h_out_base =
                        if let (Some(sc), Some(p)) = (solar_config, solar_params.as_ref()) {
                            let mut h = if sc.use_wind_speed_for_h_out {
                                sc.h_out_base_w_per_m2_k
                                    + sc.h_out_wind_coeff_w_per_m2_k_per_m_s * p.wind_speed.max(0.0)
                            } else {
                                w.h_out_w_per_m2_k
                            };
                            if sc.h_out_tilt_scale != 0.0 {
                                h *= 1.0 + sc.h_out_tilt_scale * w.normal.dz.abs();
                            }
                            h.max(1e-9)
                        } else {
                            w.h_out_w_per_m2_k.max(1e-9)
                        };
                    let h_out = super::convection::exterior_convection_h(
                        &base_config.exterior_convection_model,
                        h_out_base,
                        t_surf_ext_prev - t_out,
                        cos_tilt,
                        record.wind_speed,
                    )
                    .max(1e-9);

                    // Exterior solar + longwave
                    let mut heat_flux_net_w_per_m2 = 0.0;
                    if !w.is_ground_coupled
                        && let (Some(sc), Some(p)) = (solar_config, solar_params.as_ref())
                    {
                        if sc.include_exterior_opaque_absorption {
                            let sky_view = isotropic_sky_view_g(w.normal.dz.clamp(-1.0, 1.0));
                            let mut incident = p.diffuse_horizontal_irradiance.max(0.0) * sky_view;
                            if sun_above && p.direct_normal_irradiance > 0.0 {
                                incident += p.direct_normal_irradiance.max(0.0)
                                    * sun_dir.dot(&w.normal).max(0.0);
                            }
                            if incident > 0.0 {
                                let a = opaque_absorptance_for_path(
                                    &w.path,
                                    base_config.material_library.as_ref(),
                                    sc.default_opaque_absorptance,
                                );
                                if a > 0.0 {
                                    heat_flux_net_w_per_m2 = incident * a;
                                }
                            }
                        }
                        if sc.include_exterior_longwave_exchange {
                            let eps = sc.exterior_opaque_emissivity.clamp(0.0, 1.0);
                            if eps > 0.0 {
                                let n_z = w.normal.dz.clamp(-1.0, 1.0);
                                let sky_view = 0.5 * (1.0 + n_z);
                                let ground_view = 1.0 - sky_view;
                                let t_air_k = (p.outdoor_air_temperature_c + 273.15).max(1.0);
                                let l_sky = p.horizontal_infrared_radiation.max(0.0);
                                let eps_g = sc.ground_emissivity.clamp(0.0, 1.0);
                                let l_ground = eps_g * SIGMA_G * t_air_k.powi(4);
                                let incoming = sky_view * l_sky + ground_view * l_ground;
                                let outgoing = SIGMA_G * t_air_k.powi(4);
                                heat_flux_net_w_per_m2 += eps * (incoming - outgoing);
                            }
                        }
                    }

                    // Exterior K_eff = series(h_out*A, K_ext_face)
                    let wall_topo_idx = g_wall_conds.len();
                    let k_ext_face = topo.walls[wall_topo_idx].ext_boundary_face_k;
                    let h_out_a = h_out * w.area_m2;
                    let ext_k_eff = if k_ext_face > 0.0 {
                        1.0 / (1.0 / h_out_a + 1.0 / k_ext_face)
                    } else {
                        h_out_a
                    };

                    let ext_t_drive = if w.is_ground_coupled {
                        base_config.ground_temperature_c.unwrap_or(t_out)
                    } else {
                        t_out
                    };

                    // Exterior source (split: alpha fraction enters wall via ConvectiveWithFlux logic)
                    let alpha = if k_ext_face > 0.0 && h_out_a > 0.0 {
                        k_ext_face / (k_ext_face + h_out_a)
                    } else {
                        0.0
                    };
                    let ext_source_w = alpha * heat_flux_net_w_per_m2 * w.area_m2;

                    // Interior surface source (transmitted solar + radiant gains)
                    let mut int_source_w = 0.0;
                    if let Some(ref src) = interior_sources_walls_w_by_polygon_uid
                        && let Some(&w_src) = src.get(&w.polygon_uid)
                    {
                        int_source_w += w_src;
                    }
                    if let Some(ref src) = interior_sources_walls_w_by_zone_uid
                        && let Some(&w_total) = src.get(&w.zone_uid)
                    {
                        let a_walls: f64 = fvm_walls
                            .iter()
                            .filter(|fw| fw.area_m2 > 0.0)
                            .map(|fw| fw.area_m2)
                            .sum();
                        if a_walls > 0.0 {
                            int_source_w += w_total * (w.area_m2 / a_walls);
                        }
                    }
                    if w.normal.dz <= -0.5
                        && let Some(ref src) = floor_beam_w_by_zone_uid
                        && let Some(&beam_w) = src.get(&w.zone_uid)
                        && beam_w != 0.0
                        && floor_area > 0.0
                    {
                        int_source_w += beam_w * (w.area_m2 / floor_area);
                    }

                    g_wall_conds.push(WallStepConditions {
                        ext_k_eff,
                        ext_t_drive,
                        ext_source_w,
                        h_conv: h_in_base,
                        h_total: h_in_base,
                        int_source_w,
                    });
                }

                // Internal mass surfaces
                for m in internal_mass_surfaces.iter() {
                    let t_surf_prev = m.solver.interior_surface_temp();
                    let h_base = super::convection::interior_convection_h(
                        &base_config.interior_convection_model,
                        t_surf_prev - t_air_prev,
                        m.cos_tilt,
                        m.h_min_iso,
                    )
                    .max(1e-9);

                    // Interior surface source
                    let mut int_source_w = 0.0;
                    if let Some(ref src) = interior_sources_mass_w_by_zone_uid
                        && let Some(&w_total) = src.get(&m.zone_uid)
                    {
                        let a_mass: f64 = internal_mass_surfaces
                            .iter()
                            .filter(|s| s.area_m2 > 0.0)
                            .map(|s| s.area_m2)
                            .sum();
                        if a_mass > 0.0 {
                            int_source_w += w_total * (m.area_m2 / a_mass);
                        }
                    }

                    // For adiabatic mass: ext_k_eff=0, for TwoSided: h_conv*A coupling
                    let ext_k_eff = 0.0; // handled by topology (adiabatic or surface node)
                    g_wall_conds.push(WallStepConditions {
                        ext_k_eff,
                        ext_t_drive: t_air_prev,
                        ext_source_w: 0.0,
                        h_conv: h_base,
                        h_total: h_base,
                        int_source_w,
                    });
                }

                // Air conditions (per-zone)
                let g_air_conds: Vec<AirStepConditions> = if num_zones == 1 {
                    let direct_gains_w = gains_air_w + solar_total_air_w + q_ground;
                    vec![AirStepConditions {
                        outdoor_temp_c: t_out,
                        direct_gains_w,
                    }]
                } else {
                    // Multi-zone: distribute internal gains proportional to zone volume,
                    // solar gains are already zone-specific via glazing transmissions.
                    let total_direct = gains_air_w + solar_total_air_w + q_ground;
                    (0..num_zones)
                        .map(|zi| {
                            let frac = if volume > 0.0 {
                                zone_volumes[zi] / volume
                            } else {
                                1.0 / num_zones as f64
                            };
                            AirStepConditions {
                                outdoor_temp_c: t_out,
                                direct_gains_w: total_direct * frac,
                            }
                        })
                        .collect()
                };

                // Radiation conditions (from view factors + ScriptF)
                let g_radiation = if let Some(vf) = &view_factor_data {
                    // h_rad_base = 4σT³ (without ε — ScriptF includes emissivity)
                    let h_r_base = super::view_factors::linearized_h_rad_base(t_air_prev);

                    let n_interior_surfaces = fvm_walls.len()
                        + internal_mass_surfaces.len()
                        + global_ss_surface_map.len();
                    let mut areas = Vec::with_capacity(n_interior_surfaces);
                    for w in fvm_walls.iter() {
                        areas.push(w.area_m2);
                    }
                    for m in internal_mass_surfaces.iter() {
                        areas.push(m.area_m2);
                    }
                    for &ss_idx in &global_ss_surface_map {
                        areas.push(ss_exterior_surfaces[ss_idx].area_m2);
                    }

                    // Extract F matrix from view factor data
                    let mut f_matrix = vec![0.0; n_interior_surfaces * n_interior_surfaces];

                    // Build a handle-to-index map for our interior surface ordering
                    let mut handle_order: Vec<SurfaceHandle> = Vec::new();
                    for w in fvm_walls.iter() {
                        handle_order.push(SurfaceHandle::Polygon(w.polygon_uid.clone()));
                    }
                    for m in internal_mass_surfaces.iter() {
                        handle_order.push(SurfaceHandle::InternalMass {
                            index: m.mass_index,
                        });
                    }
                    for &ss_idx in &global_ss_surface_map {
                        handle_order.push(SurfaceHandle::Polygon(
                            ss_exterior_surfaces[ss_idx].polygon_uid.clone(),
                        ));
                    }

                    for zone_vf in &vf.zones {
                        for i in 0..zone_vf.n {
                            let handle_i = &zone_vf.surfaces[i].handle;
                            let Some(gi) = handle_order.iter().position(|h| h == handle_i) else {
                                continue;
                            };
                            for j in 0..zone_vf.n {
                                let handle_j = &zone_vf.surfaces[j].handle;
                                let Some(gj) = handle_order.iter().position(|h| h == handle_j)
                                else {
                                    continue;
                                };
                                f_matrix[gi * n_interior_surfaces + gj] =
                                    zone_vf.f_matrix[i * zone_vf.n + j];
                            }
                        }
                    }

                    // Compute ScriptF (Hottel gray enclosure exchange factors)
                    let emissivities = vec![base_config.interior_emissivity; n_interior_surfaces];
                    let script_f = super::view_factors::compute_script_f(
                        &f_matrix,
                        &emissivities,
                        n_interior_surfaces,
                    );

                    Some(RadiationConditions {
                        h_rad: h_r_base,
                        surface_areas: areas,
                        f_matrix: script_f,
                        n_surfaces: n_interior_surfaces,
                    })
                } else {
                    None
                };

                // Build wall_infos for global solve
                let g_wall_infos: Vec<FvmWallInfo> = {
                    let mut infos = Vec::new();
                    for w in fvm_walls.iter() {
                        let zi = zone_uid_to_idx.get(&w.zone_uid).copied().unwrap_or(0);
                        infos.push(FvmWallInfo {
                            solver: &w.solver,
                            zone_idx: zi,
                            area_m2: w.area_m2,
                            has_exterior_surface: false,
                            exterior_adiabatic: false,
                            exterior_zone_idx: None,
                        });
                    }
                    for m in internal_mass_surfaces.iter() {
                        let zi = zone_uid_to_idx.get(&m.zone_uid).copied().unwrap_or(0);
                        infos.push(FvmWallInfo {
                            solver: &m.solver,
                            zone_idx: zi,
                            area_m2: m.area_m2,
                            has_exterior_surface: matches!(
                                m.boundary,
                                InternalMassBoundary::TwoSided
                            ),
                            exterior_adiabatic: matches!(
                                m.boundary,
                                InternalMassBoundary::OneSidedAdiabatic
                            ),
                            exterior_zone_idx: None,
                        });
                    }
                    for (iw_idx, iw) in interzone_fvm_walls.iter().enumerate() {
                        let zi = zone_uid_to_idx.get(&iw.zone_uid).copied().unwrap_or(0);
                        let ext_zi = zone_uid_to_idx
                            .get(&interzone_exterior_zone_uids[iw_idx])
                            .copied()
                            .unwrap_or(0);
                        infos.push(FvmWallInfo {
                            solver: &iw.solver,
                            zone_idx: zi,
                            area_m2: iw.area_m2,
                            has_exterior_surface: true,
                            exterior_adiabatic: false,
                            exterior_zone_idx: Some(ext_zi),
                        });
                    }
                    infos
                };

                // Build per-zone HVAC for global solve
                let g_per_zone_hvac: Vec<HvacIdealLoads> =
                    if let Some(ref pzh) = options.per_zone_hvac {
                        pzh.clone()
                    } else {
                        vec![hour_hvac.clone(); num_zones]
                    };

                // Wall conditions for inter-zone partition walls
                for iw in interzone_fvm_walls.iter() {
                    let t_surf_prev = iw.solver.interior_surface_temp();
                    let h_in_base = super::convection::interior_convection_h(
                        &base_config.interior_convection_model,
                        t_surf_prev - t_air_prev,
                        iw.normal.dz,
                        iw.h_min_iso_interior,
                    )
                    .max(1e-9);
                    g_wall_conds.push(WallStepConditions {
                        ext_k_eff: 0.0,
                        ext_t_drive: t_air_prev,
                        ext_source_w: 0.0,
                        h_conv: h_in_base,
                        h_total: h_in_base,
                        int_source_w: 0.0,
                    });
                }

                // Solve
                let result = global_solve::step_global(
                    topo,
                    &mut global_temps,
                    &g_wall_conds,
                    &g_air_conds,
                    g_radiation.as_ref(),
                    &g_per_zone_hvac,
                    dt_s,
                    &g_wall_infos,
                );

                // Write back temperatures to wall solvers
                for (wi, w) in fvm_walls.iter_mut().enumerate() {
                    let wall_topo = &topo.walls[wi];
                    let cell_temps = &global_temps
                        [wall_topo.cell_offset..wall_topo.cell_offset + wall_topo.n_cells];
                    w.solver.set_temperatures(cell_temps);
                    // Surface temp from the surface node
                    w.cached_interior_surface_temp_c =
                        global_temps[topo.surfaces[wall_topo.surface_idx].global_idx];
                }
                let n_fvm = fvm_walls.len();
                for (mi, m) in internal_mass_surfaces.iter_mut().enumerate() {
                    let wall_topo = &topo.walls[n_fvm + mi];
                    let cell_temps = &global_temps
                        [wall_topo.cell_offset..wall_topo.cell_offset + wall_topo.n_cells];
                    m.solver.set_temperatures(cell_temps);
                    m.cached_interior_surface_temp_c =
                        global_temps[topo.surfaces[wall_topo.surface_idx].global_idx];
                }
                let n_fvm_plus_mass = n_fvm + internal_mass_surfaces.len();
                for (iw_idx, iw) in interzone_fvm_walls.iter_mut().enumerate() {
                    let wall_topo = &topo.walls[n_fvm_plus_mass + iw_idx];
                    let cell_temps = &global_temps
                        [wall_topo.cell_offset..wall_topo.cell_offset + wall_topo.n_cells];
                    iw.solver.set_temperatures(cell_temps);
                    iw.cached_interior_surface_temp_c =
                        global_temps[topo.surfaces[wall_topo.surface_idx].global_idx];
                }
                for (steady_idx, &ss_idx) in global_ss_surface_map.iter().enumerate() {
                    if let Some(ss_topo) = topo.steady_surfaces.get(steady_idx) {
                        ss_exterior_surfaces[ss_idx].cached_interior_surface_temp_c =
                            global_temps[ss_topo.global_idx];
                    }
                }

                // Update air model temperature (use zone 0 for backward compat)
                let t_air_new = global_temps[topo.air_nodes[0].global_idx];
                model.zone_temperature = t_air_new;

                // Sum heating/cooling across all zones
                let q_heating: f64 = result.heating_w_per_zone.iter().sum();
                let q_cooling: f64 = result.cooling_w_per_zone.iter().sum();
                (q_heating, q_cooling)
            };

            hour_heating_wh += heating_power * dt_h;
            hour_cooling_wh += cooling_power * dt_h;
            hour_peak_heating_w = hour_peak_heating_w.max(heating_power);
            hour_peak_cooling_w = hour_peak_cooling_w.max(cooling_power);
        }

        if report {
            hourly_heating.push(hour_heating_wh);
            hourly_cooling.push(hour_cooling_wh);

            let t_zone = model_1r1c.zone_temperature;
            hourly_zone_temp.push(t_zone);
            if t_zone < min_zone_temp {
                min_zone_temp = t_zone;
            }
            if t_zone > max_zone_temp {
                max_zone_temp = t_zone;
            }

            // Per-zone temperature tracking
            for air in global_topology.air_nodes.iter() {
                let zi = air.zone_idx;
                let t_zi = global_temps[air.global_idx];
                per_zone_hourly_temp[zi].push(t_zi);
                if t_zi < per_zone_min_temp[zi] {
                    per_zone_min_temp[zi] = t_zi;
                }
                if t_zi > per_zone_max_temp[zi] {
                    per_zone_max_temp[zi] = t_zi;
                }
            }

            annual_heating += hour_heating_wh;
            annual_cooling += hour_cooling_wh;
            peak_heating = peak_heating.max(hour_peak_heating_w);
            peak_cooling = peak_cooling.max(hour_peak_cooling_w);

            let month_idx = (record.month as usize).saturating_sub(1).min(11);
            monthly_heating[month_idx] += hour_heating_wh;
            monthly_cooling[month_idx] += hour_cooling_wh;
        }
    };

    for hour_idx in 0..warmup_hours {
        let record = &weather.records[hour_idx];
        let hour_hvac = make_hour_hvac(hvac, options, hour_idx);
        let hour_ach = options
            .hourly_infiltration_ach
            .as_ref()
            .and_then(|v| v.get(hour_idx).copied())
            .unwrap_or(base_infiltration_ach);
        simulate_hour(hour_idx, record, false, &hour_hvac, hour_ach);
    }
    for (hour_idx, record) in weather.records.iter().enumerate() {
        let hour_hvac = make_hour_hvac(hvac, options, hour_idx);
        let hour_ach = options
            .hourly_infiltration_ach
            .as_ref()
            .and_then(|v| v.get(hour_idx).copied())
            .unwrap_or(base_infiltration_ach);
        simulate_hour(hour_idx, record, true, &hour_hvac, hour_ach);
    }

    if min_zone_temp == f64::MAX {
        min_zone_temp = 0.0;
    }
    if max_zone_temp == f64::MIN {
        max_zone_temp = 0.0;
    }

    let to_kwh = 1.0 / 1000.0;

    // Fix up per-zone min/max sentinel values.
    for zi in 0..num_zones {
        if per_zone_min_temp[zi] == f64::MAX {
            per_zone_min_temp[zi] = 0.0;
        }
        if per_zone_max_temp[zi] == f64::MIN {
            per_zone_max_temp[zi] = 0.0;
        }
    }

    AnnualResult {
        hourly_heating,
        hourly_cooling,
        annual_heating_kwh: annual_heating * to_kwh,
        annual_cooling_kwh: annual_cooling * to_kwh,
        peak_heating,
        peak_cooling,
        monthly_heating_kwh: monthly_heating.map(|v| v * to_kwh),
        monthly_cooling_kwh: monthly_cooling.map(|v| v * to_kwh),
        min_zone_temp_c: min_zone_temp,
        max_zone_temp_c: max_zone_temp,
        hourly_zone_temp_c: hourly_zone_temp,
        num_zones,
        per_zone_hourly_temp_c: per_zone_hourly_temp,
        per_zone_min_temp_c: per_zone_min_temp,
        per_zone_max_temp_c: per_zone_max_temp,
    }
}

/// Builds a per-hour `HvacIdealLoads` from the base HVAC and optional schedule overrides.
fn make_hour_hvac(
    base_hvac: &HvacIdealLoads,
    options: &TransientSimulationOptions,
    hour_idx: usize,
) -> HvacIdealLoads {
    let h_set = options
        .hourly_heating_setpoint
        .as_ref()
        .and_then(|v| v.get(hour_idx).copied())
        .unwrap_or(base_hvac.heating_setpoint);
    let c_set = options
        .hourly_cooling_setpoint
        .as_ref()
        .and_then(|v| v.get(hour_idx).copied())
        .unwrap_or(base_hvac.cooling_setpoint);
    HvacIdealLoads {
        heating_setpoint: h_set,
        cooling_setpoint: c_set,
        ..*base_hvac
    }
}

/// Runs a multi-zone transient annual simulation with zone air nodes and ideal HVAC.
///
/// Differences vs [`run_transient_simulation`]:
/// - one air temperature state per `Zone` (coupled through inter-zone partitions),
/// - exterior transmission excludes internal interfaces by construction,
/// - HVAC is applied per-zone using ideal setpoints (unlimited capacity).
///
/// Gains are distributed across zones proportional to zone volume (placeholder policy
/// until per-zone schedules and solar distributions are provided by upstream modules).
pub fn run_multizone_transient_simulation(
    building: &Building,
    base_config: &ThermalConfig,
    weather: &WeatherData,
    hvac: &HvacIdealLoads,
    gains_profile: Option<&InternalGainsProfile>,
    solar_config: Option<&SolarGainConfig>,
) -> anyhow::Result<MultiZoneAnnualResult> {
    let index = SurfaceIndex::new(building);
    let boundaries = ThermalBoundaries::classify(building, &index);
    let (mut fvm_walls, fvm_skip_polygons) =
        collect_fvm_exterior_walls(building, base_config, &index, &boundaries, solar_config);
    let mut fvm_gains_by_zone_uid: std::collections::HashMap<UID, f64> =
        std::collections::HashMap::new();

    let network = if base_config.use_fvm_walls && !fvm_skip_polygons.is_empty() {
        ThermalNetwork::build_with_ignored_exterior_polygons(
            building,
            base_config,
            &index,
            &boundaries,
            &fvm_skip_polygons,
        )
    } else {
        ThermalNetwork::build(building, base_config, &index, &boundaries)
    };

    let mut model = MultiZoneAirModel::new(
        building,
        &network,
        base_config.infiltration_ach,
        base_config.thermal_capacity_j_per_m3_k,
        base_config.indoor_temperature,
    );
    let zone_uid_to_idx: std::collections::HashMap<UID, usize> = model
        .zone_uids()
        .iter()
        .enumerate()
        .map(|(i, uid)| (uid.clone(), i))
        .collect();

    let zones = building.zones();
    let zone_volumes_m3: Vec<f64> = zones.iter().map(|z| z.volume()).collect();
    let total_volume: f64 = zone_volumes_m3.iter().sum();

    let num_hours = weather.num_hours();
    let n_zones = zones.len();

    let mut hourly_zone_temperatures_c = vec![Vec::with_capacity(num_hours); n_zones];
    let mut hourly_zone_heating_w = vec![Vec::with_capacity(num_hours); n_zones];
    let mut hourly_zone_cooling_w = vec![Vec::with_capacity(num_hours); n_zones];

    let mut hourly_heating = Vec::with_capacity(num_hours);
    let mut hourly_cooling = Vec::with_capacity(num_hours);

    let mut annual_heating = 0.0;
    let mut annual_cooling = 0.0;
    let mut peak_heating = 0.0_f64;
    let mut peak_cooling = 0.0_f64;
    let mut monthly_heating = [0.0; 12];
    let mut monthly_cooling = [0.0; 12];

    for (hour_idx, record) in weather.records.iter().enumerate() {
        let gains_internal = gains_profile
            .map(|p| p.gains_at(hour_idx))
            .unwrap_or(base_config.internal_gains);
        let mut solar_params: Option<SolarHourParams> = None;
        let solar_by_zone = match solar_config {
            Some(sc) => {
                let params = SolarHourParams {
                    outdoor_air_temperature_c: record.dry_bulb_temperature,
                    global_horizontal_irradiance: record.global_horizontal_radiation,
                    direct_normal_irradiance: record.direct_normal_radiation,
                    diffuse_horizontal_irradiance: record.diffuse_horizontal_radiation,
                    horizontal_infrared_radiation: record.horizontal_infrared_radiation,
                    wind_speed: record.wind_speed,
                    day_of_year: day_of_year(record.month, record.day),
                    local_time_hours: record.hour as f64 - 0.5,
                    latitude: weather.latitude,
                    longitude: weather.longitude,
                    timezone: weather.timezone,
                };
                solar_params = Some(params);
                let mut solar = compute_solar_gains_per_zone_with_materials(
                    building,
                    &params,
                    sc,
                    base_config.material_library.as_ref(),
                );
                if sc.include_exterior_opaque_absorption || sc.include_exterior_longwave_exchange {
                    let opaque = compute_exterior_opaque_sol_air_gains_by_zone_w(
                        building,
                        base_config,
                        &index,
                        &boundaries,
                        &params,
                        sc,
                        Some(&fvm_skip_polygons),
                    );
                    for (z, q) in opaque {
                        *solar.entry(z).or_insert(0.0) += q;
                    }
                }
                solar
            }
            None => std::collections::HashMap::new(),
        };
        let mut gains_by_zone = vec![0.0; n_zones];
        if total_volume > 1e-14 {
            for i in 0..n_zones {
                let internal_i = gains_internal * (zone_volumes_m3[i] / total_volume);
                let solar_i = solar_by_zone
                    .get(&model.zone_uids()[i])
                    .cloned()
                    .unwrap_or(0.0);
                gains_by_zone[i] = internal_i + solar_i;
            }
        }

        if !fvm_walls.is_empty() {
            let temps = model.temperatures_c().to_vec();
            step_fvm_exterior_walls_fill_gains_by_zone_uid(
                &mut fvm_walls,
                base_config,
                solar_config,
                solar_params.as_ref(),
                None,
                None,
                None,
                record.dry_bulb_temperature,
                record.wind_speed,
                |zone_uid| {
                    zone_uid_to_idx
                        .get(zone_uid)
                        .and_then(|&i| temps.get(i).copied())
                        .unwrap_or(base_config.indoor_temperature)
                },
                |zone_uid| {
                    zone_uid_to_idx
                        .get(zone_uid)
                        .and_then(|&i| temps.get(i).copied())
                        .unwrap_or(base_config.indoor_temperature)
                },
                3600.0,
                &mut fvm_gains_by_zone_uid,
                None,
                0.0,
                false,
            );

            for (zone_uid, q_wall_w) in &fvm_gains_by_zone_uid {
                if let Some(&i) = zone_uid_to_idx.get(zone_uid) {
                    gains_by_zone[i] += *q_wall_w;
                }
            }
        }

        let step = model.step(record.dry_bulb_temperature, &gains_by_zone, hvac, 3600.0)?;

        let hour_heating: f64 = step.zone_heating_w.iter().sum();
        let hour_cooling: f64 = step.zone_cooling_w.iter().sum();

        hourly_heating.push(hour_heating);
        hourly_cooling.push(hour_cooling);

        annual_heating += hour_heating;
        annual_cooling += hour_cooling;
        peak_heating = peak_heating.max(hour_heating);
        peak_cooling = peak_cooling.max(hour_cooling);

        let month_idx = (record.month as usize).saturating_sub(1).min(11);
        monthly_heating[month_idx] += hour_heating;
        monthly_cooling[month_idx] += hour_cooling;

        for z in 0..n_zones {
            hourly_zone_temperatures_c[z].push(step.zone_temperatures_c[z]);
            hourly_zone_heating_w[z].push(step.zone_heating_w[z]);
            hourly_zone_cooling_w[z].push(step.zone_cooling_w[z]);
        }
    }

    let to_kwh = 1.0 / 1000.0;
    let annual = AnnualResult {
        hourly_heating,
        hourly_cooling,
        annual_heating_kwh: annual_heating * to_kwh,
        annual_cooling_kwh: annual_cooling * to_kwh,
        peak_heating,
        peak_cooling,
        monthly_heating_kwh: monthly_heating.map(|v| v * to_kwh),
        monthly_cooling_kwh: monthly_cooling.map(|v| v * to_kwh),
        min_zone_temp_c: 0.0,
        max_zone_temp_c: 0.0,
        hourly_zone_temp_c: Vec::new(),
        num_zones: 1,
        per_zone_hourly_temp_c: Vec::new(),
        per_zone_min_temp_c: Vec::new(),
        per_zone_max_temp_c: Vec::new(),
    };

    Ok(MultiZoneAnnualResult {
        zone_uids: model.zone_uids().to_vec(),
        zone_names: model.zone_names().to_vec(),
        hourly_zone_temperatures_c,
        hourly_zone_heating_w,
        hourly_zone_cooling_w,
        annual,
    })
}

/// Runs a multi-zone *steady-state* annual simulation with zone air nodes and ideal HVAC.
///
/// This uses the same multi-zone coupling model as [`run_multizone_transient_simulation`],
/// but with zero thermal capacity (instantaneous steady-state each hour). It is useful for
/// quick annual load estimates and for debugging inter-zone coupling without dynamics.
pub fn run_multizone_steady_simulation(
    building: &Building,
    base_config: &ThermalConfig,
    weather: &WeatherData,
    hvac: &HvacIdealLoads,
    gains_profile: Option<&InternalGainsProfile>,
    solar_config: Option<&SolarGainConfig>,
) -> anyhow::Result<MultiZoneAnnualResult> {
    let index = SurfaceIndex::new(building);
    let boundaries = ThermalBoundaries::classify(building, &index);
    let network = ThermalNetwork::build(building, base_config, &index, &boundaries);

    // Zero-capacity zones: steady-state per hour.
    let mut model = MultiZoneAirModel::new(
        building,
        &network,
        base_config.infiltration_ach,
        0.0,
        base_config.indoor_temperature,
    );

    let zones = building.zones();
    let zone_volumes_m3: Vec<f64> = zones.iter().map(|z| z.volume()).collect();
    let total_volume: f64 = zone_volumes_m3.iter().sum();

    let num_hours = weather.num_hours();
    let n_zones = zones.len();

    let mut hourly_zone_temperatures_c = vec![Vec::with_capacity(num_hours); n_zones];
    let mut hourly_zone_heating_w = vec![Vec::with_capacity(num_hours); n_zones];
    let mut hourly_zone_cooling_w = vec![Vec::with_capacity(num_hours); n_zones];

    let mut hourly_heating = Vec::with_capacity(num_hours);
    let mut hourly_cooling = Vec::with_capacity(num_hours);

    let mut annual_heating = 0.0;
    let mut annual_cooling = 0.0;
    let mut peak_heating = 0.0_f64;
    let mut peak_cooling = 0.0_f64;
    let mut monthly_heating = [0.0; 12];
    let mut monthly_cooling = [0.0; 12];

    for (hour_idx, record) in weather.records.iter().enumerate() {
        let gains_internal = gains_profile.map(|p| p.gains_at(hour_idx)).unwrap_or(0.0);
        let solar_by_zone = match solar_config {
            Some(sc) => {
                let params = SolarHourParams {
                    outdoor_air_temperature_c: record.dry_bulb_temperature,
                    global_horizontal_irradiance: record.global_horizontal_radiation,
                    direct_normal_irradiance: record.direct_normal_radiation,
                    diffuse_horizontal_irradiance: record.diffuse_horizontal_radiation,
                    horizontal_infrared_radiation: record.horizontal_infrared_radiation,
                    wind_speed: record.wind_speed,
                    day_of_year: day_of_year(record.month, record.day),
                    local_time_hours: record.hour as f64 - 0.5,
                    latitude: weather.latitude,
                    longitude: weather.longitude,
                    timezone: weather.timezone,
                };
                compute_solar_gains_per_zone_with_materials(
                    building,
                    &params,
                    sc,
                    base_config.material_library.as_ref(),
                )
            }
            None => std::collections::HashMap::new(),
        };

        let mut gains_by_zone = vec![0.0; n_zones];
        if total_volume > 1e-14 {
            for i in 0..n_zones {
                let internal_i = gains_internal * (zone_volumes_m3[i] / total_volume);
                let solar_i = solar_by_zone
                    .get(&model.zone_uids()[i])
                    .cloned()
                    .unwrap_or(0.0);
                gains_by_zone[i] = internal_i + solar_i;
            }
        }

        let step = model.step(record.dry_bulb_temperature, &gains_by_zone, hvac, 3600.0)?;

        let hour_heating: f64 = step.zone_heating_w.iter().sum();
        let hour_cooling: f64 = step.zone_cooling_w.iter().sum();

        hourly_heating.push(hour_heating);
        hourly_cooling.push(hour_cooling);

        annual_heating += hour_heating;
        annual_cooling += hour_cooling;
        peak_heating = peak_heating.max(hour_heating);
        peak_cooling = peak_cooling.max(hour_cooling);

        let month_idx = (record.month as usize).saturating_sub(1).min(11);
        monthly_heating[month_idx] += hour_heating;
        monthly_cooling[month_idx] += hour_cooling;

        for z in 0..n_zones {
            hourly_zone_temperatures_c[z].push(step.zone_temperatures_c[z]);
            hourly_zone_heating_w[z].push(step.zone_heating_w[z]);
            hourly_zone_cooling_w[z].push(step.zone_cooling_w[z]);
        }
    }

    let to_kwh = 1.0 / 1000.0;
    let annual = AnnualResult {
        hourly_heating,
        hourly_cooling,
        annual_heating_kwh: annual_heating * to_kwh,
        annual_cooling_kwh: annual_cooling * to_kwh,
        peak_heating,
        peak_cooling,
        monthly_heating_kwh: monthly_heating.map(|v| v * to_kwh),
        monthly_cooling_kwh: monthly_cooling.map(|v| v * to_kwh),
        min_zone_temp_c: 0.0,
        max_zone_temp_c: 0.0,
        hourly_zone_temp_c: Vec::new(),
        num_zones: 1,
        per_zone_hourly_temp_c: Vec::new(),
        per_zone_min_temp_c: Vec::new(),
        per_zone_max_temp_c: Vec::new(),
    };

    Ok(MultiZoneAnnualResult {
        zone_uids: model.zone_uids().to_vec(),
        zone_names: model.zone_names().to_vec(),
        hourly_zone_temperatures_c,
        hourly_zone_heating_w,
        hourly_zone_cooling_w,
        annual,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sim::energy::weather::HourlyRecord;
    use crate::{Solid, Zone};

    #[test]
    fn test_annual_simulation_basic() {
        let s = Solid::from_box(5.0, 5.0, 3.0, None, "room").unwrap();
        let zone = Zone::new("z", vec![s]).unwrap();
        let building = Building::new("b", vec![zone]).unwrap();

        let config = ThermalConfig::new();
        let weather = WeatherData::synthetic("Test", 52.0, 13.0, 10.0, 12.0);

        let result = run_annual_simulation(&building, &config, &weather, None, None);

        assert_eq!(result.hourly_heating.len(), 8760);
        assert_eq!(result.hourly_cooling.len(), 8760);
        assert!(
            result.annual_heating_kwh > 0.0,
            "Should need heating in cold climate"
        );
        assert!(result.peak_heating > 0.0, "Should have peak heating load");
    }

    #[test]
    fn test_annual_with_gains() {
        let s = Solid::from_box(10.0, 10.0, 3.0, None, "room").unwrap();
        let zone = Zone::new("z", vec![s]).unwrap();
        let building = Building::new("b", vec![zone]).unwrap();

        let config = ThermalConfig::new();
        let weather = WeatherData::synthetic("Test", 52.0, 13.0, 10.0, 12.0);
        let gains = InternalGainsProfile::office(100.0);

        let result_no_gains = run_annual_simulation(&building, &config, &weather, None, None);
        let result_gains = run_annual_simulation(&building, &config, &weather, Some(&gains), None);

        // Internal gains should reduce heating demand
        assert!(
            result_gains.annual_heating_kwh < result_no_gains.annual_heating_kwh,
            "Internal gains should reduce heating: {} vs {}",
            result_gains.annual_heating_kwh,
            result_no_gains.annual_heating_kwh
        );
    }

    #[test]
    fn test_monthly_sums() {
        let s = Solid::from_box(3.0, 3.0, 3.0, None, "room").unwrap();
        let zone = Zone::new("z", vec![s]).unwrap();
        let building = Building::new("b", vec![zone]).unwrap();

        let config = ThermalConfig::new();
        let weather = WeatherData::synthetic("Test", 52.0, 13.0, 10.0, 12.0);

        let result = run_annual_simulation(&building, &config, &weather, None, None);

        let monthly_sum: f64 = result.monthly_heating_kwh.iter().sum();
        assert!(
            (monthly_sum - result.annual_heating_kwh).abs() < 1.0,
            "Monthly heating sum should equal annual: {monthly_sum} vs {}",
            result.annual_heating_kwh
        );
    }

    #[test]
    fn test_day_of_year_basic() {
        assert_eq!(day_of_year(1, 1), 1);
        assert_eq!(day_of_year(3, 1), 60); // non-leap-year convention (Jan31+Feb28+1)
        assert_eq!(day_of_year(12, 31), 365);
    }

    #[test]
    fn test_annual_simulation_with_solar_config_small_weather() {
        let s = Solid::from_box(5.0, 5.0, 3.0, None, "room").unwrap();
        let zone = Zone::new("z", vec![s]).unwrap();
        let building = Building::new("b", vec![zone]).unwrap();

        let config = ThermalConfig::new();
        let weather = WeatherData {
            location: "Mini".to_string(),
            latitude: 0.0,
            longitude: 0.0,
            timezone: 0.0,
            elevation: 0.0,
            records: vec![
                HourlyRecord {
                    month: 3,
                    day: 20,
                    hour: 12,
                    dry_bulb_temperature: 0.0,
                    relative_humidity: 50.0,
                    horizontal_infrared_radiation: 300.0,
                    global_horizontal_radiation: 0.0,
                    direct_normal_radiation: 800.0,
                    diffuse_horizontal_radiation: 100.0,
                    wind_speed: 0.0,
                    wind_direction: 0.0,
                },
                HourlyRecord {
                    month: 3,
                    day: 20,
                    hour: 1,
                    dry_bulb_temperature: 0.0,
                    relative_humidity: 50.0,
                    horizontal_infrared_radiation: 300.0,
                    global_horizontal_radiation: 0.0,
                    direct_normal_radiation: 0.0,
                    diffuse_horizontal_radiation: 50.0,
                    wind_speed: 0.0,
                    wind_direction: 0.0,
                },
            ],
        };

        let mut solar = SolarGainConfig::new();
        solar.default_shgc = 1.0;
        solar.glazing_patterns = vec!["floor".to_string()];

        let no_solar = run_annual_simulation(&building, &config, &weather, None, None);
        let with_solar = run_annual_simulation(&building, &config, &weather, None, Some(&solar));

        assert_eq!(no_solar.hourly_heating.len(), 2);
        assert_eq!(with_solar.hourly_heating.len(), 2);
        assert!(
            with_solar.hourly_heating[0] <= no_solar.hourly_heating[0],
            "Solar gains should not increase heating demand"
        );
    }

    #[test]
    fn test_multizone_steady_simulation_hits_solar_per_zone_branch() {
        let s = Solid::from_box(5.0, 5.0, 3.0, None, "room").unwrap();
        let zone = Zone::new("z", vec![s]).unwrap();
        let building = Building::new("b", vec![zone]).unwrap();

        let mut config = ThermalConfig::new();
        config.default_u_value = 1.0;
        config.infiltration_ach = 0.0;
        config.indoor_temperature = 20.0;
        config.outdoor_temperature = 20.0;

        let weather = WeatherData {
            location: "Mini".to_string(),
            latitude: 0.0,
            longitude: 0.0,
            timezone: 0.0,
            elevation: 0.0,
            records: vec![HourlyRecord {
                month: 6,
                day: 1,
                hour: 12,
                dry_bulb_temperature: 20.0,
                relative_humidity: 50.0,
                horizontal_infrared_radiation: 300.0,
                global_horizontal_radiation: 0.0,
                direct_normal_radiation: 0.0,
                diffuse_horizontal_radiation: 200.0,
                wind_speed: 0.0,
                wind_direction: 0.0,
            }],
        };

        let hvac = HvacIdealLoads::new();

        let no_solar =
            run_multizone_steady_simulation(&building, &config, &weather, &hvac, None, None)
                .unwrap();

        let mut solar = SolarGainConfig::new();
        solar.default_shgc = 1.0;
        solar.glazing_patterns = vec!["floor".to_string()];

        let with_solar = run_multizone_steady_simulation(
            &building,
            &config,
            &weather,
            &hvac,
            None,
            Some(&solar),
        )
        .unwrap();

        assert_eq!(no_solar.annual.hourly_heating.len(), 1);
        assert_eq!(with_solar.annual.hourly_heating.len(), 1);
        assert_eq!(no_solar.annual.hourly_cooling.len(), 1);
        assert_eq!(with_solar.annual.hourly_cooling.len(), 1);
        assert!(
            with_solar.annual.hourly_cooling[0] > no_solar.annual.hourly_cooling[0],
            "Solar gains should increase cooling demand in steady-state"
        );
    }

    #[test]
    fn test_transient_simulation() {
        let s = Solid::from_box(5.0, 5.0, 3.0, None, "room").unwrap();
        let zone = Zone::new("z", vec![s]).unwrap();
        let building = Building::new("b", vec![zone]).unwrap();

        let config = ThermalConfig::new();
        let weather = WeatherData::synthetic("Test", 52.0, 13.0, 10.0, 12.0);
        let hvac = HvacIdealLoads::new();

        let result = run_transient_simulation(&building, &config, &weather, &hvac, None, None);

        assert_eq!(result.hourly_heating.len(), 8760);
        assert!(
            result.annual_heating_kwh > 0.0,
            "Should need heating in cold climate"
        );
        assert!(result.peak_heating > 0.0, "Should have peak heating load");

        // Transient should have different results from steady-state due to thermal mass
        let steady = run_annual_simulation(&building, &config, &weather, None, None);
        // They won't be identical, but both should indicate heating is needed
        assert!(steady.annual_heating_kwh > 0.0);
    }

    #[test]
    fn test_transient_simulation_dt_zero_fallback_ua_and_solar() {
        let s = Solid::from_box(2.0, 2.0, 2.0, None, "room").unwrap();
        let zone = Zone::new("z", vec![s]).unwrap();
        let building = Building::new("b", vec![zone]).unwrap();

        let mut config = ThermalConfig::new();
        config.indoor_temperature = 20.0;
        config.outdoor_temperature = 20.0; // dt=0 -> triggers UA/infiltration fallback branches

        let weather = WeatherData {
            location: "Mini".to_string(),
            latitude: 0.0,
            longitude: 0.0,
            timezone: 0.0,
            elevation: 0.0,
            records: vec![
                HourlyRecord {
                    month: 6,
                    day: 1,
                    hour: 12,
                    dry_bulb_temperature: 20.0,
                    relative_humidity: 50.0,
                    horizontal_infrared_radiation: 300.0,
                    global_horizontal_radiation: 0.0,
                    direct_normal_radiation: 500.0,
                    diffuse_horizontal_radiation: 50.0,
                    wind_speed: 0.0,
                    wind_direction: 0.0,
                },
                HourlyRecord {
                    month: 6,
                    day: 1,
                    hour: 1,
                    dry_bulb_temperature: 20.0,
                    relative_humidity: 50.0,
                    horizontal_infrared_radiation: 300.0,
                    global_horizontal_radiation: 0.0,
                    direct_normal_radiation: 0.0,
                    diffuse_horizontal_radiation: 10.0,
                    wind_speed: 0.0,
                    wind_direction: 0.0,
                },
            ],
        };

        let hvac = HvacIdealLoads::new();
        let gains = InternalGainsProfile::office(100.0);
        let mut solar = SolarGainConfig::new();
        solar.default_shgc = 1.0;
        solar.glazing_patterns = vec!["floor".to_string()];

        let result = run_transient_simulation(
            &building,
            &config,
            &weather,
            &hvac,
            Some(&gains),
            Some(&solar),
        );
        assert_eq!(result.hourly_heating.len(), 2);
        assert_eq!(result.hourly_cooling.len(), 2);
    }

    #[test]
    fn test_multizone_transient_simulation_basic() {
        let s0 = Solid::from_box(2.0, 2.0, 2.0, None, "s0").unwrap();
        let s1 = Solid::from_box(2.0, 2.0, 2.0, Some((2.0, 0.0, 0.0)), "s1").unwrap();
        let z0 = Zone::new("z0", vec![s0]).unwrap();
        let z1 = Zone::new("z1", vec![s1]).unwrap();
        let building = Building::new("b", vec![z0, z1]).unwrap();

        let config = ThermalConfig::new();
        let weather = WeatherData::synthetic("Test", 52.0, 13.0, 10.0, 12.0);
        let hvac = HvacIdealLoads::new();

        let result =
            run_multizone_transient_simulation(&building, &config, &weather, &hvac, None, None)
                .unwrap();

        assert_eq!(result.zone_names.len(), 2);
        assert_eq!(result.hourly_zone_temperatures_c.len(), 2);
        assert_eq!(result.hourly_zone_temperatures_c[0].len(), 8760);
        assert_eq!(result.annual.hourly_heating.len(), 8760);
    }

    #[test]
    fn test_multizone_transient_with_solar_distribution() {
        let s0 = Solid::from_box(2.0, 2.0, 2.0, None, "s0").unwrap();
        let s1 = Solid::from_box(2.0, 2.0, 2.0, Some((2.0, 0.0, 0.0)), "s1").unwrap();
        let z0 = Zone::new("z0", vec![s0]).unwrap();
        let z1 = Zone::new("z1", vec![s1]).unwrap();
        let building = Building::new("b", vec![z0, z1]).unwrap();

        let config = ThermalConfig::new();
        let weather = WeatherData {
            location: "Mini".to_string(),
            latitude: 0.0,
            longitude: 0.0,
            timezone: 0.0,
            elevation: 0.0,
            records: vec![
                HourlyRecord {
                    month: 6,
                    day: 1,
                    hour: 12,
                    dry_bulb_temperature: 0.0,
                    relative_humidity: 50.0,
                    global_horizontal_radiation: 0.0,
                    direct_normal_radiation: 800.0,
                    diffuse_horizontal_radiation: 100.0,
                    horizontal_infrared_radiation: 300.0,
                    wind_speed: 0.0,
                    wind_direction: 0.0,
                },
                HourlyRecord {
                    month: 6,
                    day: 1,
                    hour: 1,
                    dry_bulb_temperature: 0.0,
                    relative_humidity: 50.0,
                    global_horizontal_radiation: 0.0,
                    direct_normal_radiation: 0.0,
                    diffuse_horizontal_radiation: 20.0,
                    horizontal_infrared_radiation: 300.0,
                    wind_speed: 0.0,
                    wind_direction: 0.0,
                },
            ],
        };
        let hvac = HvacIdealLoads::new();
        let gains = InternalGainsProfile::office(100.0);
        let mut solar = SolarGainConfig::new();
        solar.default_shgc = 1.0;
        solar.glazing_patterns = vec!["floor".to_string()];

        let result = run_multizone_transient_simulation(
            &building,
            &config,
            &weather,
            &hvac,
            Some(&gains),
            Some(&solar),
        )
        .unwrap();

        assert_eq!(result.zone_names.len(), 2);
        assert_eq!(result.hourly_zone_temperatures_c.len(), 2);
        assert_eq!(result.hourly_zone_temperatures_c[0].len(), 2);
        assert_eq!(result.annual.hourly_heating.len(), 2);
        assert_eq!(result.annual.hourly_cooling.len(), 2);
    }

    #[test]
    fn test_multizone_steady_simulation_basic() {
        let s0 = Solid::from_box(2.0, 2.0, 2.0, None, "s0").unwrap();
        let s1 = Solid::from_box(2.0, 2.0, 2.0, Some((2.0, 0.0, 0.0)), "s1").unwrap();
        let z0 = Zone::new("z0", vec![s0]).unwrap();
        let z1 = Zone::new("z1", vec![s1]).unwrap();
        let building = Building::new("b", vec![z0, z1]).unwrap();

        let config = ThermalConfig::new();
        let weather = WeatherData::synthetic("Test", 52.0, 13.0, 10.0, 12.0);
        let hvac = HvacIdealLoads::new();

        let result =
            run_multizone_steady_simulation(&building, &config, &weather, &hvac, None, None)
                .unwrap();

        assert_eq!(result.zone_names.len(), 2);
        assert_eq!(result.hourly_zone_temperatures_c.len(), 2);
        assert_eq!(result.hourly_zone_temperatures_c[0].len(), 8760);
        assert_eq!(result.annual.hourly_heating.len(), 8760);
    }

    // ── Physics verification tests ──────────────────────────────────────

    #[test]
    fn test_annual_energy_conservation() {
        // With a 1-hour timestep, summing hourly demands (W) over all hours is numerically
        // equivalent to the annual energy in Wh, i.e. annual_kwh * 1000.
        let s = Solid::from_box(4.0, 4.0, 3.0, None, "room").unwrap();
        let zone = Zone::new("z", vec![s]).unwrap();
        let building = Building::new("b", vec![zone]).unwrap();

        let config = ThermalConfig::new();
        let weather = WeatherData::synthetic("Test", 50.0, 10.0, 8.0, 10.0);

        let result = run_annual_simulation(&building, &config, &weather, None, None);

        let sum_heating_wh: f64 = result.hourly_heating.iter().sum();
        let sum_cooling_wh: f64 = result.hourly_cooling.iter().sum();

        assert!(
            (sum_heating_wh / 1000.0 - result.annual_heating_kwh).abs() < 0.01,
            "Hourly heating sum={:.1} Wh vs annual={:.1} kWh",
            sum_heating_wh,
            result.annual_heating_kwh
        );
        assert!(
            (sum_cooling_wh / 1000.0 - result.annual_cooling_kwh).abs() < 0.01,
            "Hourly cooling sum={:.1} Wh vs annual={:.1} kWh",
            sum_cooling_wh,
            result.annual_cooling_kwh
        );
    }

    #[test]
    fn test_peak_is_max_of_hourly() {
        let s = Solid::from_box(5.0, 5.0, 3.0, None, "room").unwrap();
        let zone = Zone::new("z", vec![s]).unwrap();
        let building = Building::new("b", vec![zone]).unwrap();

        let config = ThermalConfig::new();
        let weather = WeatherData::synthetic("Test", 52.0, 13.0, 10.0, 15.0);

        let result = run_annual_simulation(&building, &config, &weather, None, None);

        let max_heating = result
            .hourly_heating
            .iter()
            .cloned()
            .fold(0.0_f64, f64::max);
        let max_cooling = result
            .hourly_cooling
            .iter()
            .cloned()
            .fold(0.0_f64, f64::max);

        assert!(
            (result.peak_heating - max_heating).abs() < 1e-10,
            "peak_heating={} vs max(hourly)={}",
            result.peak_heating,
            max_heating
        );
        assert!(
            (result.peak_cooling - max_cooling).abs() < 1e-10,
            "peak_cooling={} vs max(hourly)={}",
            result.peak_cooling,
            max_cooling
        );
    }

    #[test]
    fn test_annual_hour_matches_steady_state_heat_balance() {
        let s = Solid::from_box(5.0, 5.0, 3.0, None, "room").unwrap();
        let zone = Zone::new("z", vec![s]).unwrap();
        let building = Building::new("b", vec![zone]).unwrap();

        let mut config = ThermalConfig::new();
        config.indoor_temperature = 20.0;

        // Synthetic weather with 0 annual amplitude (still has daily variation)
        let weather = WeatherData::synthetic("Const", 50.0, 10.0, 5.0, 0.0);

        let result = run_annual_simulation(&building, &config, &weather, None, None);

        // Annual simulation runs a steady-state heat balance every hour. Spot-check one hour.
        let record = &weather.records[0];
        let mut hour_config = config.clone();
        hour_config.outdoor_temperature = record.dry_bulb_temperature;
        hour_config.internal_gains = 0.0;
        hour_config.solar_gains = 0.0;
        let hour_result = calculate_heat_balance(&building, &hour_config);

        assert!(
            (result.hourly_heating[0] - hour_result.heating_demand).abs() < 1e-12,
            "Hour 0 heating mismatch: annual={} vs heat_balance={}",
            result.hourly_heating[0],
            hour_result.heating_demand
        );
        assert!(
            (result.hourly_cooling[0] - hour_result.cooling_demand).abs() < 1e-12,
            "Hour 0 cooling mismatch: annual={} vs heat_balance={}",
            result.hourly_cooling[0],
            hour_result.cooling_demand
        );
    }

    #[test]
    fn test_synthetic_weather_amplitude_zero_repeats_daily_cycle() {
        // Synthetic weather always has a ±3°C daily swing; with `temp_amplitude=0.0`,
        // the annual component is removed, so the day-to-day pattern repeats.
        let s = Solid::from_box(5.0, 5.0, 3.0, None, "room").unwrap();
        let zone = Zone::new("z", vec![s]).unwrap();
        let building = Building::new("b", vec![zone]).unwrap();

        let mut config = ThermalConfig::new();
        config.indoor_temperature = 20.0;

        let weather = WeatherData::synthetic("Const", 50.0, 10.0, 5.0, 0.0);
        let result = run_annual_simulation(&building, &config, &weather, None, None);

        // Compare day 1 vs day 2 (hours 0..24 vs 24..48).
        for h in 0..24 {
            let day1 = result.hourly_heating[h];
            let day2 = result.hourly_heating[h + 24];
            assert!(
                (day1 - day2).abs() < 1e-6,
                "Hour-of-day {h}: day1={day1}, day2={day2}"
            );
        }

        // All heating values should be non-negative
        for (i, &val) in result.hourly_heating.iter().enumerate() {
            assert!(
                val >= 0.0,
                "Hour {i}: heating should be non-negative, got {val}"
            );
        }
    }

    #[test]
    fn test_transient_energy_conservation() {
        // For the transient simulation, hourly sums should also match annual kWh.
        let s = Solid::from_box(5.0, 5.0, 3.0, None, "room").unwrap();
        let zone = Zone::new("z", vec![s]).unwrap();
        let building = Building::new("b", vec![zone]).unwrap();

        let config = ThermalConfig::new();
        let weather = WeatherData::synthetic("Test", 52.0, 13.0, 10.0, 12.0);
        let hvac = HvacIdealLoads::new();

        let result = run_transient_simulation(&building, &config, &weather, &hvac, None, None);

        let sum_heating: f64 = result.hourly_heating.iter().sum();
        let sum_cooling: f64 = result.hourly_cooling.iter().sum();

        assert!(
            (sum_heating / 1000.0 - result.annual_heating_kwh).abs() < 0.01,
            "Transient hourly heating sum mismatch"
        );
        assert!(
            (sum_cooling / 1000.0 - result.annual_cooling_kwh).abs() < 0.01,
            "Transient hourly cooling sum mismatch"
        );
    }

    #[test]
    fn test_more_insulation_less_heating() {
        // Better insulated building should need less heating.
        let s = Solid::from_box(5.0, 5.0, 3.0, None, "room").unwrap();
        let zone = Zone::new("z", vec![s]).unwrap();
        let building = Building::new("b", vec![zone]).unwrap();

        let weather = WeatherData::synthetic("Cold", 55.0, 10.0, 5.0, 12.0);

        let mut config_poor = ThermalConfig::new();
        config_poor.default_u_value = 3.0;

        let mut config_good = ThermalConfig::new();
        config_good.default_u_value = 0.5;

        let result_poor = run_annual_simulation(&building, &config_poor, &weather, None, None);
        let result_good = run_annual_simulation(&building, &config_good, &weather, None, None);

        assert!(
            result_good.annual_heating_kwh < result_poor.annual_heating_kwh,
            "Better insulation should reduce heating: good={}, poor={}",
            result_good.annual_heating_kwh,
            result_poor.annual_heating_kwh
        );
        // With U=0.5 vs U=3.0, the ratio of transmission losses is 6:1.
        // Infiltration is unchanged, so the total ratio is less than 6,
        // but the well-insulated building should use significantly less.
        assert!(
            result_good.annual_heating_kwh < result_poor.annual_heating_kwh * 0.5,
            "U=0.5 should use well under half the energy of U=3.0"
        );
    }

    #[test]
    fn test_opaque_solar_gains() {
        let s = Solid::from_box(5.0, 5.0, 3.0, None, "room").unwrap();
        let zone = Zone::new("z", vec![s]).unwrap();
        let building = Building::new("b", vec![zone]).unwrap();

        let weather = WeatherData {
            location: "Solar".to_string(),
            latitude: 30.0,
            longitude: 0.0,
            timezone: 0.0,
            elevation: 0.0,
            records: vec![HourlyRecord {
                month: 6,
                day: 21,
                hour: 12,
                dry_bulb_temperature: 0.0,
                relative_humidity: 50.0,
                horizontal_infrared_radiation: 300.0,
                global_horizontal_radiation: 0.0,
                direct_normal_radiation: 800.0,
                diffuse_horizontal_radiation: 200.0,
                wind_speed: 0.0,
                wind_direction: 0.0,
            }],
        };
        let hvac = HvacIdealLoads::new();

        // Without opaque solar
        let config = ThermalConfig::new();
        let result_no_opaque =
            run_transient_simulation(&building, &config, &weather, &hvac, None, None);

        // With opaque solar absorption
        let mut solar = SolarGainConfig::new();
        solar.include_exterior_opaque_absorption = true;
        solar.default_opaque_absorptance = 0.9;
        let result_opaque =
            run_transient_simulation(&building, &config, &weather, &hvac, None, Some(&solar));

        // Opaque solar absorption should reduce heating demand
        // (exterior solar heats the surfaces which conducts inward)
        assert!(
            result_opaque.hourly_heating[0] <= result_no_opaque.hourly_heating[0],
            "Opaque solar should reduce heating: with={}, without={}",
            result_opaque.hourly_heating[0],
            result_no_opaque.hourly_heating[0]
        );
    }

    #[test]
    fn test_transient_with_warmup_hours() {
        let s = Solid::from_box(5.0, 5.0, 3.0, None, "room").unwrap();
        let zone = Zone::new("z", vec![s]).unwrap();
        let building = Building::new("b", vec![zone]).unwrap();

        let weather = WeatherData::synthetic("Test", 52.0, 13.0, 10.0, 12.0);
        let hvac = HvacIdealLoads::new();

        let config = ThermalConfig::new();

        let opts_warmup = TransientSimulationOptions {
            warmup_hours: 48,
            substeps_per_hour: 1,
            ..Default::default()
        };

        let result_warmup = run_transient_simulation_with_options(
            &building,
            &config,
            &weather,
            &hvac,
            None,
            None,
            &opts_warmup,
        );

        // Should produce 8760 hours of results (warmup doesn't add to reported)
        assert_eq!(result_warmup.hourly_heating.len(), 8760);
        assert!(result_warmup.annual_heating_kwh > 0.0);
    }

    #[test]
    fn test_transient_with_solar_and_opaque_absorption() {
        let s = Solid::from_box(5.0, 5.0, 3.0, None, "room").unwrap();
        let zone = Zone::new("z", vec![s]).unwrap();
        let building = Building::new("b", vec![zone]).unwrap();

        let weather = WeatherData {
            location: "Solar".to_string(),
            latitude: 30.0,
            longitude: 0.0,
            timezone: 0.0,
            elevation: 0.0,
            records: vec![
                HourlyRecord {
                    month: 6,
                    day: 21,
                    hour: 12,
                    dry_bulb_temperature: 10.0,
                    relative_humidity: 50.0,
                    horizontal_infrared_radiation: 350.0,
                    global_horizontal_radiation: 0.0,
                    direct_normal_radiation: 800.0,
                    diffuse_horizontal_radiation: 200.0,
                    wind_speed: 3.0,
                    wind_direction: 0.0,
                },
                HourlyRecord {
                    month: 6,
                    day: 21,
                    hour: 1,
                    dry_bulb_temperature: 5.0,
                    relative_humidity: 50.0,
                    horizontal_infrared_radiation: 300.0,
                    global_horizontal_radiation: 0.0,
                    direct_normal_radiation: 0.0,
                    diffuse_horizontal_radiation: 10.0,
                    wind_speed: 1.0,
                    wind_direction: 0.0,
                },
            ],
        };
        let hvac = HvacIdealLoads::new();
        let config = ThermalConfig::new();

        // Enable all sol-air features
        let mut solar = SolarGainConfig::new();
        solar.include_exterior_opaque_absorption = true;
        solar.include_exterior_longwave_exchange = true;
        solar.use_wind_speed_for_h_out = true;

        let result =
            run_transient_simulation(&building, &config, &weather, &hvac, None, Some(&solar));
        assert_eq!(result.hourly_heating.len(), 2);
    }

    #[test]
    fn test_fvm_h_out_tilt_scale_applies_without_wind_model() {
        use crate::sim::energy::construction::WallConstruction;
        use crate::sim::energy::convection::ExteriorConvectionModel;
        use crate::sim::materials::Layer;

        let construction = WallConstruction::single_layer(
            "test",
            Layer {
                name: "layer".to_string(),
                thickness: 0.2,
                conductivity: 1.0,
                density: 2000.0,
                specific_heat: 900.0,
            },
        );
        let mesh = build_1d_mesh(&construction, 10.0);
        let solver = FvmWallSolver::new(mesh, 20.0);
        let wall = FvmExteriorWall {
            zone_uid: UID::from("zone"),
            polygon_uid: UID::new(),
            path: "zone/solid/roof/poly".to_string(),
            area_m2: 10.0,
            normal: crate::Vector::new(0.0, 0.0, 1.0),
            is_ground_coupled: false,
            h_out_w_per_m2_k: 8.0,
            h_min_iso_interior: 8.0,
            cached_interior_surface_temp_c: solver.interior_surface_temp(),
            solver,
        };

        let mut walls_no_tilt = vec![wall.clone()];
        let mut walls_with_tilt = vec![wall];

        let mut config = ThermalConfig::new();
        config.exterior_convection_model = ExteriorConvectionModel::Fixed;
        let params = SolarHourParams {
            outdoor_air_temperature_c: 0.0,
            global_horizontal_irradiance: 0.0,
            direct_normal_irradiance: 0.0,
            diffuse_horizontal_irradiance: 0.0,
            horizontal_infrared_radiation: 0.0,
            wind_speed: 0.0,
            day_of_year: 172,
            local_time_hours: 12.0,
            latitude: 0.0,
            longitude: 0.0,
            timezone: 0.0,
        };

        let mut solar_no_tilt = SolarGainConfig::new();
        solar_no_tilt.use_wind_speed_for_h_out = false;
        solar_no_tilt.h_out_tilt_scale = 0.0;
        let mut solar_with_tilt = solar_no_tilt.clone();
        solar_with_tilt.h_out_tilt_scale = 1.0;

        let mut gains_no_tilt = std::collections::HashMap::new();
        let mut gains_with_tilt = std::collections::HashMap::new();

        step_fvm_exterior_walls_fill_gains_by_zone_uid(
            &mut walls_no_tilt,
            &config,
            Some(&solar_no_tilt),
            Some(&params),
            None,
            None,
            None,
            0.0,
            0.0,
            |_| 20.0,
            |_| 20.0,
            3600.0,
            &mut gains_no_tilt,
            None,
            0.0,
            false,
        );
        step_fvm_exterior_walls_fill_gains_by_zone_uid(
            &mut walls_with_tilt,
            &config,
            Some(&solar_with_tilt),
            Some(&params),
            None,
            None,
            None,
            0.0,
            0.0,
            |_| 20.0,
            |_| 20.0,
            3600.0,
            &mut gains_with_tilt,
            None,
            0.0,
            false,
        );

        let q_no_tilt = *gains_no_tilt.get(&UID::from("zone")).unwrap_or(&0.0);
        let q_with_tilt = *gains_with_tilt.get(&UID::from("zone")).unwrap_or(&0.0);

        assert!(
            q_no_tilt < 0.0,
            "Expected heat loss to zone air, got {q_no_tilt}"
        );
        assert!(
            q_with_tilt.abs() > q_no_tilt.abs(),
            "Tilt scaling should increase |q| with fixed h_out; no_tilt={q_no_tilt}, with_tilt={q_with_tilt}"
        );
    }

    #[test]
    fn test_fvm_view_factor_air_gain_is_convective_only() {
        use crate::sim::energy::construction::WallConstruction;
        use crate::sim::energy::convection::InteriorConvectionModel;
        use crate::sim::materials::Layer;

        let construction = WallConstruction::single_layer(
            "test",
            Layer {
                name: "layer".to_string(),
                thickness: 0.2,
                conductivity: 1.0,
                density: 2000.0,
                specific_heat: 900.0,
            },
        );
        let area_m2 = 10.0;
        let mesh = build_1d_mesh(&construction, area_m2);
        let solver = FvmWallSolver::new(mesh, 20.0);
        let polygon_uid = UID::new();
        let mut walls = vec![FvmExteriorWall {
            zone_uid: UID::from("zone"),
            polygon_uid: polygon_uid.clone(),
            path: "zone/solid/wall/poly".to_string(),
            area_m2,
            normal: crate::Vector::new(0.0, 0.0, 1.0),
            is_ground_coupled: false,
            h_out_w_per_m2_k: 8.0,
            h_min_iso_interior: 8.0,
            cached_interior_surface_temp_c: solver.interior_surface_temp(),
            solver,
        }];

        let mut config = ThermalConfig::new();
        config.interior_convection_model = InteriorConvectionModel::Fixed(3.0);
        config.interior_emissivity = 0.9;

        let mut gains = std::collections::HashMap::new();
        let mut mrt_map = std::collections::HashMap::new();
        let t_air = 20.0;
        let t_mrt = 100.0;
        let h_r_uniform = 10.0;
        mrt_map.insert(SurfaceHandle::Polygon(polygon_uid), t_mrt);

        step_fvm_exterior_walls_fill_gains_by_zone_uid(
            &mut walls,
            &config,
            None,
            None,
            None,
            None,
            None,
            20.0,
            0.0,
            |_| t_air,
            |_| t_air,
            3600.0,
            &mut gains,
            Some(&mrt_map),
            h_r_uniform,
            false,
        );

        let q_to_air = *gains.get(&UID::from("zone")).unwrap_or(&0.0);
        let t_surf = walls[0].cached_interior_surface_temp_c;
        let h_conv = 3.0;
        let q_conv_expected = h_conv * (t_surf - t_air) * area_m2;
        let tol = 1e-9 * q_conv_expected.abs().max(1.0);
        assert!(
            (q_to_air - q_conv_expected).abs() <= tol,
            "air gain should be convection-only in VF mode; got {q_to_air}, expected {q_conv_expected}"
        );

        // Guard against regression to leaky total-coupling accounting.
        let h_rad = config.interior_emissivity * h_r_uniform;
        let h_total = h_conv + h_rad;
        let t_eff = (h_conv * t_air + h_rad * t_mrt) / h_total;
        let q_leaky = h_total * (t_surf - t_eff) * area_m2;
        assert!(
            (q_to_air - q_leaky).abs() > 1e-3,
            "air gain must not include radiative exchange term; leaky formula matched unexpectedly"
        );
    }

    #[test]
    fn test_multizone_transient_with_solar_per_zone_distribution() {
        // Exercise compute_exterior_opaque_sol_air_gains_by_zone_w
        let s0 = Solid::from_box(3.0, 3.0, 3.0, None, "s0").unwrap();
        let s1 = Solid::from_box(3.0, 3.0, 3.0, Some((3.0, 0.0, 0.0)), "s1").unwrap();
        let z0 = Zone::new("z0", vec![s0]).unwrap();
        let z1 = Zone::new("z1", vec![s1]).unwrap();
        let building = Building::new("b", vec![z0, z1]).unwrap();

        let weather = WeatherData {
            location: "S".to_string(),
            latitude: 30.0,
            longitude: 0.0,
            timezone: 0.0,
            elevation: 0.0,
            records: vec![HourlyRecord {
                month: 6,
                day: 21,
                hour: 12,
                dry_bulb_temperature: 10.0,
                relative_humidity: 50.0,
                horizontal_infrared_radiation: 350.0,
                global_horizontal_radiation: 0.0,
                direct_normal_radiation: 800.0,
                diffuse_horizontal_radiation: 200.0,
                wind_speed: 2.0,
                wind_direction: 0.0,
            }],
        };
        let hvac = HvacIdealLoads::new();
        let config = ThermalConfig::new();

        let mut solar = SolarGainConfig::new();
        solar.include_exterior_opaque_absorption = true;
        solar.include_exterior_longwave_exchange = true;
        solar.use_wind_speed_for_h_out = true;

        let result = run_multizone_transient_simulation(
            &building,
            &config,
            &weather,
            &hvac,
            None,
            Some(&solar),
        )
        .unwrap();

        assert_eq!(result.zone_names.len(), 2);
        assert_eq!(result.hourly_zone_temperatures_c[0].len(), 1);
    }

    #[test]
    fn test_opaque_absorptance_for_path_function() {
        use crate::sim::materials::{Material, MaterialLibrary, OpticalMaterial};

        // No library
        let a = opaque_absorptance_for_path("z/s/w/p", None, 0.7);
        assert!((a - 0.7).abs() < 1e-12);

        // Library with no matching material
        let mut lib = MaterialLibrary::new();
        lib.add(Material::new("mat"));
        lib.assign("wall", "mat");
        let a = opaque_absorptance_for_path("unmatched", Some(&lib), 0.7);
        assert!((a - 0.7).abs() < 1e-12);

        // Library with matching material but no optical
        let a = opaque_absorptance_for_path("z/s/wall/p", Some(&lib), 0.7);
        assert!((a - 0.7).abs() < 1e-12);

        // Library with optical properties
        let mut mat = Material::new("dark");
        mat.optical = Some(OpticalMaterial {
            name: "dark".to_string(),
            diffuse_reflectance: [0.2, 0.2, 0.2],
            specular_reflectance: [0.0; 3],
            transmittance: [0.0; 3],
        });
        lib.add(mat);
        lib.assign("roof", "dark");
        let a = opaque_absorptance_for_path("z/s/roof/p", Some(&lib), 0.7);
        assert!((a - 0.8).abs() < 1e-12, "Expected 0.8, got {a}");
    }
}
