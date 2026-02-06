use anyhow::Result;
use rerun as rr;

use crate::sim::acoustics::impulse_response::ImpulseResponse;
use crate::sim::acoustics::receiver::Receiver;
use crate::sim::materials::NUM_OCTAVE_BANDS;

const SESSION_NAME: &str = "Building3d";

/// Draws receivers as colored spheres in the 3D scene.
pub fn draw_receivers(session: &rr::RecordingStream, receivers: &[Receiver]) -> Result<()> {
    let mut pts = Vec::new();
    let mut radii = Vec::new();
    let mut colors = Vec::new();

    for receiver in receivers {
        pts.push(receiver.position);
        radii.push(receiver.radius as f32);
        colors.push(rr::Color(rr::Rgba32::from_linear_unmultiplied_rgba_f32(
            0.0, 0.5, 1.0, 0.6,
        )));
    }

    session.log_static(
        format!("{}/receivers", SESSION_NAME),
        &rr::Points3D::new(pts).with_radii(radii).with_colors(colors),
    )?;

    Ok(())
}

/// Draws an impulse response decay curve as scalar time series.
///
/// Logs the broadband Schroeder decay curve and per-band energy levels.
pub fn draw_impulse_response(
    session: &rr::RecordingStream,
    ir: &ImpulseResponse,
    receiver_name: &str,
) -> Result<()> {
    let decay = ir.schroeder_decay_broadband();

    // Log broadband Schroeder decay
    for (i, &level) in decay.iter().enumerate() {
        let time_ms = (i as f64 * ir.time_resolution * 1000.0) as i64;
        session.set_time_sequence("time_ms", time_ms);

        let db = if level > 1e-30 {
            10.0 * level.log10()
        } else {
            -60.0
        };

        session.log(
            format!("{}/ir/{}/schroeder_db", SESSION_NAME, receiver_name),
            &rr::Scalars::new([db]),
        )?;
    }

    // Log per-band energy
    let band_labels = ["125Hz", "250Hz", "500Hz", "1kHz", "2kHz", "4kHz"];
    for band in 0..NUM_OCTAVE_BANDS.min(ir.bands.len()) {
        let band_decay = ir.schroeder_decay(band);
        for (i, &level) in band_decay.iter().enumerate() {
            let time_ms = (i as f64 * ir.time_resolution * 1000.0) as i64;
            session.set_time_sequence("time_ms", time_ms);

            let db = if level > 1e-30 {
                10.0 * level.log10()
            } else {
                -60.0
            };

            let label = if band < band_labels.len() {
                band_labels[band]
            } else {
                "band"
            };
            session.log(
                format!("{}/ir/{}/{}", SESSION_NAME, receiver_name, label),
                &rr::Scalars::new([db]),
            )?;
        }
    }

    Ok(())
}
