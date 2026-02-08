//! BRAS CR2 validation benchmark.
//!
//! Loads the CR2 (small seminar room) geometry from BRAS, assigns fitted
//! material absorption/scattering data, runs frequency-dependent acoustic
//! ray tracing, and compares computed room acoustic parameters (EDT, T20,
//! C80, D50) against in-situ measurements.
//!
//! Reference: BRAS (Benchmark for Room Acoustical Simulation), CC BY-SA 4.0
//! https://depositonce.tu-berlin.de/items/version/64

use anyhow::{Context, Result};
use std::path::Path;

use building3d::Point;
use building3d::io::{Ac3dCoordSystem, read_ac3d};
use building3d::sim::acoustics::impulse_response::ImpulseResponse;
use building3d::sim::acoustics::metrics::RoomAcousticReport;
use building3d::sim::acoustics::receiver::Receiver;
use building3d::sim::materials::{
    AcousticMaterial, Material, MaterialLibrary, NUM_OCTAVE_BANDS, OCTAVE_BAND_FREQUENCIES,
};
use building3d::sim::rays::{AcousticMode, Simulation, SimulationConfig};

/// Third-octave band indices (in the 31-band 20 Hz–20 kHz CSV) for the
/// 6 octave bands used by building3d (125, 250, 500, 1000, 2000, 4000 Hz).
const OCTAVE_INDICES: [usize; 6] = [8, 11, 14, 17, 20, 23];

/// Parse a BRAS material CSV file.
///
/// Format: 3 lines of 31 comma-separated values each.
///   Line 1: frequencies (ignored)
///   Line 2: absorption coefficients
///   Line 3: scattering coefficients
///
/// Returns the 6 octave-band values extracted at `OCTAVE_INDICES`.
fn parse_material_csv(path: &Path) -> Result<([f64; 6], [f64; 6])> {
    let content = std::fs::read_to_string(path)
        .with_context(|| format!("Failed to read material CSV: {}", path.display()))?;
    let lines: Vec<&str> = content.lines().collect();
    if lines.len() < 3 {
        anyhow::bail!(
            "Material CSV must have at least 3 lines: {}",
            path.display()
        );
    }

    let parse_line = |line: &str| -> Result<[f64; 6]> {
        let values: Vec<f64> = line
            .split(',')
            .map(|s| s.trim().parse::<f64>())
            .collect::<std::result::Result<Vec<_>, _>>()
            .with_context(|| "Failed to parse CSV values")?;
        if values.len() < 31 {
            anyhow::bail!("Expected 31 values, got {}", values.len());
        }
        let mut out = [0.0; 6];
        for (i, &idx) in OCTAVE_INDICES.iter().enumerate() {
            out[i] = values[idx];
        }
        Ok(out)
    };

    let absorption = parse_line(lines[1])?;
    let scattering = parse_line(lines[2])?;
    Ok((absorption, scattering))
}

fn main() -> Result<()> {
    let base = Path::new("validation/bras");
    let ac3d_path = base.join("informed_sim/RavenModels/scene9.ac");
    let csv_dir = base.join("surface_descriptions/3 Surface descriptions/_csv/fitted_estimates");

    // ── 1. Load geometry ──────────────────────────────────────────────
    println!("Loading CR2 geometry from {}...", ac3d_path.display());
    let (_mat_names, building) = read_ac3d(&ac3d_path, "seminar_room", Ac3dCoordSystem::ZUp)?;
    println!("  Loaded building with {} zones", building.zones().len());

    // ── 2. Parse material CSVs → MaterialLibrary ──────────────────────
    let material_names = ["concrete", "windows", "ceiling", "plaster", "floor"];
    let mut lib = MaterialLibrary::new();
    for name in &material_names {
        let csv_path = csv_dir.join(format!("mat_CR2_{}.csv", name));
        let (abs, scat) = parse_material_csv(&csv_path)
            .with_context(|| format!("Failed to parse material '{}'", name))?;
        println!(
            "  Material {:>10}: abs=[{:.3}, {:.3}, {:.3}, {:.3}, {:.3}, {:.3}]",
            name, abs[0], abs[1], abs[2], abs[3], abs[4], abs[5]
        );
        lib.add(Material::new(name).with_acoustic(AcousticMaterial::new(name, abs, scat)));
        // Wall names in the AC3D model match material names, so substring
        // match on the material name will find the correct wall.
        lib.assign(name, name);
    }

    // ── 3. Source and receiver positions ───────────────────────────────
    // RAVEN Y-up: (x, y_up, z) → building3d Z-up: (x, -z, y_up)
    let sources = [
        ("LS1", Point::new(0.931, -2.547, 0.723)),
        ("LS2", Point::new(0.119, 2.880, 0.723)),
    ];
    let receivers = [
        ("MP1", Point::new(-0.993, 1.426, 1.230)),
        ("MP2", Point::new(0.439, -0.147, 1.230)),
        ("MP3", Point::new(1.361, -0.603, 1.230)),
        ("MP4", Point::new(-1.110, -0.256, 1.230)),
        ("MP5", Point::new(-0.998, -1.409, 1.230)),
    ];

    let absorber_positions: Vec<Point> = receivers.iter().map(|(_, p)| *p).collect();

    // ── 4. Simulation parameters ──────────────────────────────────────
    let time_step = 2.5e-5; // 25 μs
    let num_rays = 50_000;
    let max_time = 3.0; // 3 seconds (RAVEN uses 2.8 s)
    let num_steps = (max_time / time_step) as usize;
    let absorber_radius = 0.8;
    // Receiver time resolution for histograms (RAVEN uses 2 ms)
    let receiver_dt = 0.002;

    println!("\nSimulation parameters:");
    println!("  Rays: {}", num_rays);
    println!(
        "  Steps: {} (dt={:.1e} s, total={:.1} s)",
        num_steps, time_step, max_time
    );
    println!("  Absorber radius: {} m", absorber_radius);
    println!("  Receiver time resolution: {} ms", receiver_dt * 1000.0);

    // ── 5. Run simulation for each source ─────────────────────────────
    // Accumulate per-pair reports for averaging
    let mut all_reports: Vec<(String, RoomAcousticReport)> = Vec::new();

    for (src_name, src_pos) in sources.iter() {
        println!("\n{}", "=".repeat(60));
        println!(
            "Source {} at ({:.3}, {:.3}, {:.3})",
            src_name, src_pos.x, src_pos.y, src_pos.z
        );
        println!("{}", "=".repeat(60));

        let mut config = SimulationConfig::new();
        config.acoustic_mode = AcousticMode::FrequencyDependent;
        config.material_library = Some(lib.clone());
        config.enable_air_absorption = true;
        config.search_transparent = false;
        config.source = *src_pos;
        config.absorbers = absorber_positions.clone();
        config.absorber_radius = absorber_radius;
        config.num_rays = num_rays;
        config.num_steps = num_steps;
        config.time_step = time_step;
        config.store_ray_history = false;
        config.min_alive_fraction = 0.01; // Stop when <1% of rays are alive

        println!("  Running simulation...");
        let sim = Simulation::new(&building, config)?;
        let result = sim.run();
        println!("  Done. {} steps computed.", result.hits.len());

        // Extract metrics for each receiver
        let band_hits = result
            .band_hits
            .as_ref()
            .expect("FrequencyDependent mode should produce band_hits");

        for (ri, (recv_name, recv_pos)) in receivers.iter().enumerate() {
            let mut receiver = Receiver::new(*recv_pos, absorber_radius, receiver_dt, max_time);

            // Feed band hits from simulation result
            for (step, absorber_hits) in band_hits.iter().enumerate() {
                let time = step as f64 * time_step;
                if ri < absorber_hits.len() {
                    let energy: f64 = absorber_hits[ri].iter().sum();
                    if energy > 0.0 {
                        receiver.record_band_hit(time, &absorber_hits[ri]);
                    }
                }
            }

            let ir = ImpulseResponse::from_receiver(&receiver);
            let report = RoomAcousticReport::from_ir(&ir);

            let pair_name = format!("{}->{}", src_name, recv_name);
            println!(
                "\n  {} (total energy: {:.6}):",
                pair_name,
                receiver.total_energy()
            );
            print_report(&report);

            all_reports.push((pair_name, report));
        }
    }

    // ── 6. Average across all pairs and compare with measurements ─────
    println!("\n\n{}", "=".repeat(80));
    println!("COMPARISON: Simulated (averaged) vs Measured");
    println!("{}", "=".repeat(80));

    // Measured reference data (from SUMMARY.md, averaged across all source-receiver pairs)
    #[rustfmt::skip]
    let measured_edt: [f64; 6]  = [1.501, 1.306, 2.018, 1.919, 1.758, 1.649];
    #[rustfmt::skip]
    let measured_t20: [f64; 6]  = [1.453, 1.345, 2.018, 1.934, 1.715, 1.604];
    #[rustfmt::skip]
    let measured_c80: [f64; 6]  = [0.551, 0.902, -1.862, -0.606, -0.483, 0.164];
    #[rustfmt::skip]
    let measured_d50: [f64; 6]  = [0.295, 0.394, 0.278, 0.321, 0.356, 0.366]; // converted from %

    // Average simulated metrics
    let _n = all_reports.len() as f64;
    let mut avg_edt = [0.0f64; NUM_OCTAVE_BANDS];
    let mut avg_t20 = [0.0f64; NUM_OCTAVE_BANDS];
    let mut avg_c80 = [0.0f64; NUM_OCTAVE_BANDS];
    let mut avg_d50 = [0.0f64; NUM_OCTAVE_BANDS];
    let mut cnt_edt = [0usize; NUM_OCTAVE_BANDS];
    let mut cnt_t20 = [0usize; NUM_OCTAVE_BANDS];
    let mut cnt_c80 = [0usize; NUM_OCTAVE_BANDS];
    let mut cnt_d50 = [0usize; NUM_OCTAVE_BANDS];

    for (_, report) in &all_reports {
        for b in 0..NUM_OCTAVE_BANDS {
            if let Some(v) = report.edt[b] {
                avg_edt[b] += v;
                cnt_edt[b] += 1;
            }
            if let Some(v) = report.rt60[b] {
                // rt60 uses T30 with T20 fallback; we want T20 specifically
                // but RoomAcousticReport stores rt60. Use metrics::t20 directly
                // if needed. For now use rt60 as proxy.
                avg_t20[b] += v;
                cnt_t20[b] += 1;
            }
            if let Some(v) = report.c80[b] {
                avg_c80[b] += v;
                cnt_c80[b] += 1;
            }
            if let Some(v) = report.d50[b] {
                avg_d50[b] += v;
                cnt_d50[b] += 1;
            }
        }
    }

    println!(
        "\n  {:>8} | {:>8} {:>8} {:>6} | {:>8} {:>8} {:>6} | {:>8} {:>8} {:>6} | {:>7} {:>7} {:>6}",
        "Freq",
        "Sim EDT",
        "Meas",
        "Err%",
        "Sim T20",
        "Meas",
        "Err%",
        "Sim C80",
        "Meas",
        "Err",
        "Sim D50",
        "Meas",
        "Err"
    );
    println!("  {}", "-".repeat(120));

    for b in 0..NUM_OCTAVE_BANDS {
        let freq = OCTAVE_BAND_FREQUENCIES[b];
        let sim_edt = if cnt_edt[b] > 0 {
            avg_edt[b] / cnt_edt[b] as f64
        } else {
            f64::NAN
        };
        let sim_t20 = if cnt_t20[b] > 0 {
            avg_t20[b] / cnt_t20[b] as f64
        } else {
            f64::NAN
        };
        let sim_c80 = if cnt_c80[b] > 0 {
            avg_c80[b] / cnt_c80[b] as f64
        } else {
            f64::NAN
        };
        let sim_d50 = if cnt_d50[b] > 0 {
            avg_d50[b] / cnt_d50[b] as f64
        } else {
            f64::NAN
        };

        let edt_err = if measured_edt[b].abs() > 0.01 {
            (sim_edt - measured_edt[b]) / measured_edt[b] * 100.0
        } else {
            f64::NAN
        };
        let t20_err = if measured_t20[b].abs() > 0.01 {
            (sim_t20 - measured_t20[b]) / measured_t20[b] * 100.0
        } else {
            f64::NAN
        };
        let c80_err = sim_c80 - measured_c80[b]; // absolute dB error
        let d50_err = (sim_d50 - measured_d50[b]) * 100.0; // absolute percentage-point error

        println!(
            "  {:>6.0} Hz | {:>8.3} {:>8.3} {:>+5.1}% | {:>8.3} {:>8.3} {:>+5.1}% | {:>7.2} dB {:>7.2} {:>+5.1} | {:>6.1}% {:>6.1}% {:>+5.1}",
            freq,
            sim_edt,
            measured_edt[b],
            edt_err,
            sim_t20,
            measured_t20[b],
            t20_err,
            sim_c80,
            measured_c80[b],
            c80_err,
            sim_d50 * 100.0,
            measured_d50[b] * 100.0,
            d50_err,
        );
    }

    println!("\n  Note: T20 column shows RT60 (T30 with T20 fallback) from simulation.");
    println!("  {} source-receiver pairs averaged.", all_reports.len());

    Ok(())
}

fn print_report(report: &RoomAcousticReport) {
    println!(
        "    {:>8} {:>8} {:>8} {:>9} {:>7}",
        "Freq", "EDT(s)", "RT60(s)", "C80(dB)", "D50"
    );
    for b in 0..NUM_OCTAVE_BANDS {
        let freq = report.frequencies[b];
        let edt_s = report.edt[b].map_or("--".to_string(), |v| format!("{:.3}", v));
        let rt60_s = report.rt60[b].map_or("--".to_string(), |v| format!("{:.3}", v));
        let c80_s = report.c80[b].map_or("--".to_string(), |v| format!("{:.2}", v));
        let d50_s = report.d50[b].map_or("--".to_string(), |v| format!("{:.1}%", v * 100.0));
        println!(
            "    {:>6.0} Hz {:>8} {:>8} {:>9} {:>7}",
            freq, edt_s, rt60_s, c80_s, d50_s
        );
    }
}
