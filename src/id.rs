use uuid::Uuid;

/// Returns a random UUID
pub fn random_id() -> String {
    Uuid::new_v4().to_string()
}
