//! Unique identifier type used throughout the building hierarchy.

use serde::{Deserialize, Serialize};
use uuid::Uuid;

/// A thin wrapper around a UUID v4 string, used as a stable identity for every
/// entity in the building hierarchy.
#[derive(Eq, PartialEq, Hash, Debug, Clone, Serialize, Deserialize)]
pub struct UID(String);

impl From<&str> for UID {
    fn from(value: &str) -> Self {
        Self(value.to_string())
    }
}

impl From<String> for UID {
    fn from(value: String) -> Self {
        Self(value)
    }
}

impl Default for UID {
    fn default() -> Self {
        Self::new()
    }
}

impl UID {
    /// Generates a fresh random UUID v4 identifier.
    pub fn new() -> Self {
        Self(Self::random())
    }

    /// Returns the underlying UUID string slice.
    pub fn as_str(&self) -> &str {
        &self.0
    }

    fn random() -> String {
        Uuid::new_v4().to_string()
    }
}
