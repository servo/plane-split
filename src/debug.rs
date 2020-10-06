use crate::{Polygon, Splitter};
use serde::{Serialize, Serializer};

/// Serialized work for a plane splitter.
pub struct Dump<T, U, A> {
    /// input polygons
    input: Vec<Polygon<T, U, A>>,
    /// view used to sort
    view: euclid::Vector3D<T, U>,
    /// split polygons
    output: Vec<Polygon<T, U, A>>,
}

impl<T: Serialize, U, A: Serialize> Serialize for Dump<T, U, A> {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        use serde::ser::SerializeStruct;
        let mut me = serializer.serialize_struct("Dump", 3)?;
        me.serialize_field("input", &self.input)?;
        me.serialize_field("view", &self.view)?;
        me.serialize_field("output", &self.output)?;
        me.end()
    }
}

/// Debug layer that records the interface into a dump.
pub struct DebugLayer<T, U, A, Z> {
    /// Actual plane splitting implementation.
    inner: Z,
    /// Dump of the work.
    dump: Dump<T, U, A>,
}

impl<T: Default, U, A, Z> DebugLayer<T, U, A, Z> {
    /// Create a new debug layer.
    pub fn new(inner: Z) -> Self {
        DebugLayer {
            inner,
            dump: Dump {
                input: Vec::new(),
                view: Default::default(),
                output: Vec::new(),
            },
        }
    }

    /// Get the current work dump.
    pub fn dump(&self) -> &Dump<T, U, A> {
        &self.dump
    }
}

impl<T: Clone, U, A: Copy, Z: Splitter<T, U, A>> Splitter<T, U, A> for DebugLayer<T, U, A, Z> {
    fn reset(&mut self) {
        self.dump.input.clear();
        self.inner.reset();
    }

    /// Add a new polygon and return a slice of the subdivisions
    /// that avoid collision with any of the previously added polygons.
    fn add(&mut self, polygon: Polygon<T, U, A>) {
        self.dump.input.push(polygon.clone());
        self.inner.add(polygon);
    }

    /// Sort the produced polygon set by the ascending distance across
    /// the specified view vector. Return the sorted slice.
    fn sort(&mut self, view: euclid::Vector3D<T, U>) -> &[Polygon<T, U, A>] {
        self.dump.view = view.clone();
        let sorted = self.inner.sort(view);
        self.dump.output.clear();
        self.dump.output.extend_from_slice(sorted);
        sorted
    }
}
