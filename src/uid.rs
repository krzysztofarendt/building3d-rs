use serde::{Deserialize, Serialize};
use uuid::Uuid;

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
    pub fn new() -> Self {
        Self(Self::random())
    }

    pub fn as_str(&self) -> &str {
        &self.0
    }

    fn random() -> String {
        Uuid::new_v4().to_string()
    }
}
