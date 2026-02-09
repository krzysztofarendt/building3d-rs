use crate::Building;

use crate::sim::index::SurfaceIndex;

/// Shared read-only context passed to simulation modules.
///
/// Domain-specific metadata (e.g. "exterior wall", "inter-zone partition",
/// "glazing") should not be stored on geometry types. Instead, modules derive
/// the metadata they need from the geometry and from user-provided assignments,
/// and store it in their own overlay data structures (often keyed by `UID`).
pub struct SimContext<'a> {
    pub building: &'a Building,
    pub surface_index: &'a SurfaceIndex,
}

impl<'a> SimContext<'a> {
    pub fn new(building: &'a Building, surface_index: &'a SurfaceIndex) -> Self {
        Self {
            building,
            surface_index,
        }
    }
}
