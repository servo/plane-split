extern crate euclid;
extern crate plane_split;

use euclid::TypedPoint3D;
use plane_split::{NaiveSplitter, Splitter, _make_grid};

#[test]
fn naive() {
    let count = 2;
    let polys = _make_grid(count);
    let mut splitter = NaiveSplitter::new();
    let result = splitter.solve(&polys, TypedPoint3D::new(0.0, 0.0, 1.0));
    assert_eq!(result.len(), count + count*count + count*count*count);
}
