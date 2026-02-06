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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_str() {
        let uid = UID::from("test-uid");
        assert_eq!(uid.as_str(), "test-uid");
    }

    #[test]
    fn test_from_string() {
        let uid = UID::from(String::from("test-uid-2"));
        assert_eq!(uid.as_str(), "test-uid-2");
    }

    #[test]
    fn test_default() {
        let uid = UID::default();
        // Default generates a UUID, should be non-empty
        assert!(!uid.as_str().is_empty());
    }

    #[test]
    fn test_as_str() {
        let uid = UID::new();
        let s = uid.as_str();
        assert!(!s.is_empty());
        // UUID v4 format: 8-4-4-4-12 hex chars
        assert_eq!(s.len(), 36);
    }

    #[test]
    fn test_uniqueness() {
        let uid1 = UID::new();
        let uid2 = UID::new();
        assert_ne!(uid1, uid2);
    }
}
