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

/// Minimum dB floor for Schroeder decay visualization.
const DB_FLOOR: f64 = -60.0;

/// Draws an impulse response decay curve as scalar time series.
///
/// Logs the broadband Schroeder decay curve and per-band energy levels.
/// Values are clamped to -60 dB floor to keep the chart readable.
pub fn draw_impulse_response(
    session: &rr::RecordingStream,
    ir: &ImpulseResponse,
    receiver_name: &str,
) -> Result<()> {
    let decay = ir.schroeder_decay_broadband();

    // Log broadband Schroeder decay on the "step" timeline so it appears
    // alongside the ray animation in Rerun.
    // Note: schroeder_decay_broadband() already returns values in dB.
    for (i, &db) in decay.iter().enumerate() {
        let db = db.max(DB_FLOOR);
        session.set_time_sequence("step", i as i64);
        session.log(
            format!("{}/ir/{}/schroeder_db", SESSION_NAME, receiver_name),
            &rr::Scalars::new([db]),
        )?;
    }

    // Log per-band Schroeder decay (already in dB)
    let band_labels = ["125Hz", "250Hz", "500Hz", "1kHz", "2kHz", "4kHz"];
    for band in 0..NUM_OCTAVE_BANDS.min(ir.bands.len()) {
        let band_decay = ir.schroeder_decay(band);
        for (i, &db) in band_decay.iter().enumerate() {
            let db = db.max(DB_FLOOR);
            session.set_time_sequence("step", i as i64);
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Point;

    fn buffered_session() -> rr::RecordingStream {
        rr::RecordingStreamBuilder::new("test").buffered().unwrap()
    }

    #[test]
    fn test_draw_receivers_is_ok() {
        let session = buffered_session();
        let receivers = vec![Receiver::new(Point::new(0.0, 0.0, 0.0), 0.25, 0.001, 0.1)];
        draw_receivers(&session, &receivers).unwrap();
    }

    #[test]
    fn test_draw_impulse_response_is_ok() {
        let session = buffered_session();

        let mut receiver = Receiver::new(Point::new(0.0, 0.0, 0.0), 0.25, 0.001, 0.2);
        // Add some energy early so Schroeder integration has non-zero values.
        receiver.record_scalar_hit(0.001, 1.0);
        receiver.record_scalar_hit(0.010, 0.2);
        let ir = ImpulseResponse::from_receiver(&receiver);

        draw_impulse_response(&session, &ir, "r0").unwrap();
    }
}
