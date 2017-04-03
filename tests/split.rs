extern crate plane_split;

use plane_split::{NaiveSplitter, Splitter, _make_grid};

#[test]
fn naive() {
    let count = 2;
    let polys = _make_grid(count);
    let mut splitter = NaiveSplitter::new();
    let result = splitter.solve(&polys);
    assert_eq!(result.len(), count + count*count + count*count*count);
}
