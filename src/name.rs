/// Types that expose a comparable name.
pub trait HasName {
    fn get_name(&self) -> &str;
}

// Delegate HasName to references (and smart pointers if useful)
impl<T: HasName + ?Sized> HasName for &T {
    fn get_name(&self) -> &str { (*self).get_name() }
}
impl<T: HasName + ?Sized> HasName for Box<T> {
    fn get_name(&self) -> &str { (**self).get_name() }
}
impl<T: HasName + ?Sized> HasName for std::rc::Rc<T> {
    fn get_name(&self) -> &str { (**self).get_name() }
}
impl<T: HasName + ?Sized> HasName for std::sync::Arc<T> {
    fn get_name(&self) -> &str { (**self).get_name() }
}

/// Sorting helpers for slices of `T: HasName`.
pub trait SortByName {
    /// Stable, ascending sort by `name()`.
    fn sort_by_name(&mut self);
}

impl<T: HasName> SortByName for [T] {
    fn sort_by_name(&mut self) {
        // `sort_by` is stable since Rust 1.2; compares &str by Unicode scalar values.
        self.sort_by(|a, b| a.get_name().cmp(b.get_name()));
    }
}
