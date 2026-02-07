/// Types that expose a comparable name.
pub trait HasName {
    fn get_name(&self) -> &str;
}

// Delegate HasName to references (and smart pointers if useful)
impl<T: HasName + ?Sized> HasName for &T {
    fn get_name(&self) -> &str {
        (*self).get_name()
    }
}
impl<T: HasName + ?Sized> HasName for Box<T> {
    fn get_name(&self) -> &str {
        (**self).get_name()
    }
}
impl<T: HasName + ?Sized> HasName for std::rc::Rc<T> {
    fn get_name(&self) -> &str {
        (**self).get_name()
    }
}
impl<T: HasName + ?Sized> HasName for std::sync::Arc<T> {
    fn get_name(&self) -> &str {
        (**self).get_name()
    }
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

#[cfg(test)]
mod tests {
    use super::*;
    use std::rc::Rc;
    use std::sync::Arc;

    struct Named(String);
    impl HasName for Named {
        fn get_name(&self) -> &str {
            &self.0
        }
    }

    #[test]
    fn test_has_name_box() {
        let item: Box<Named> = Box::new(Named("hello".to_string()));
        assert_eq!(item.get_name(), "hello");
    }

    #[test]
    fn test_has_name_rc() {
        let item: Rc<Named> = Rc::new(Named("world".to_string()));
        assert_eq!(item.get_name(), "world");
    }

    #[test]
    fn test_has_name_arc() {
        let item: Arc<Named> = Arc::new(Named("arc_name".to_string()));
        assert_eq!(item.get_name(), "arc_name");
    }

    #[test]
    fn test_sort_by_name() {
        let mut items = vec![
            Named("charlie".to_string()),
            Named("alice".to_string()),
            Named("bob".to_string()),
        ];
        items.as_mut_slice().sort_by_name();
        assert_eq!(items[0].get_name(), "alice");
        assert_eq!(items[1].get_name(), "bob");
        assert_eq!(items[2].get_name(), "charlie");
    }
}
