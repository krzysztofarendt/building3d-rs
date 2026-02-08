use std::any::{Any, TypeId};
use std::collections::HashMap;

/// Typed message/value storage used to connect independent simulation modules.
///
/// Values are keyed by their concrete type, so personalized simulations can
/// introduce new message types without modifying a central enum.
#[derive(Default)]
pub struct Bus {
    values: HashMap<TypeId, Box<dyn Any>>,
}

impl Bus {
    pub fn new() -> Self {
        Self::default()
    }

    /// Inserts or replaces the stored value of type `T`.
    pub fn put<T: 'static>(&mut self, value: T) {
        self.values.insert(TypeId::of::<T>(), Box::new(value));
    }

    /// Inserts the value only if no value of type `T` is already present.
    ///
    /// Returns an error if a value of that type already exists on the Bus.
    /// Use this to enforce single-producer contracts in composed pipelines.
    pub fn put_unique<T: 'static>(&mut self, value: T) -> anyhow::Result<()> {
        use std::collections::hash_map::Entry;
        match self.values.entry(TypeId::of::<T>()) {
            Entry::Vacant(e) => {
                e.insert(Box::new(value));
                Ok(())
            }
            Entry::Occupied(_) => {
                anyhow::bail!(
                    "Bus already contains a value of type `{}`; \
                     a composed pipeline must have exactly one producer per type",
                    std::any::type_name::<T>(),
                );
            }
        }
    }

    /// Gets a reference to the stored value of type `T`, if present.
    pub fn get<T: 'static>(&self) -> Option<&T> {
        self.values
            .get(&TypeId::of::<T>())
            .and_then(|v| v.downcast_ref::<T>())
    }

    /// Gets a mutable reference to the stored value of type `T`, if present.
    pub fn get_mut<T: 'static>(&mut self) -> Option<&mut T> {
        self.values
            .get_mut(&TypeId::of::<T>())
            .and_then(|v| v.downcast_mut::<T>())
    }

    /// Removes and returns the stored value of type `T`, if present.
    pub fn take<T: 'static>(&mut self) -> Option<T> {
        self.values
            .remove(&TypeId::of::<T>())
            .and_then(|v| v.downcast::<T>().ok())
            .map(|b| *b)
    }
}
