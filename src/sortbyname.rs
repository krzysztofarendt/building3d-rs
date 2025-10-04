/// Types that expose a comparable name.
pub trait HasName {
    fn name(&self) -> &str;
}

// Delegate HasName to references (and smart pointers if useful)
impl<T: HasName + ?Sized> HasName for &T {
    fn name(&self) -> &str { (*self).name() }
}
impl<T: HasName + ?Sized> HasName for Box<T> {
    fn name(&self) -> &str { (**self).name() }
}
impl<T: HasName + ?Sized> HasName for std::rc::Rc<T> {
    fn name(&self) -> &str { (**self).name() }
}
impl<T: HasName + ?Sized> HasName for std::sync::Arc<T> {
    fn name(&self) -> &str { (**self).name() }
}

/// Sorting helpers for slices of `T: HasName`.
pub trait SortByName {
    /// Stable, ascending sort by `name()`.
    fn sort_by_name(&mut self);
}

impl<T: HasName> SortByName for [T] {
    fn sort_by_name(&mut self) {
        // `sort_by` is stable since Rust 1.2; compares &str by Unicode scalar values.
        self.sort_by(|a, b| a.name().cmp(b.name()));
    }
}
