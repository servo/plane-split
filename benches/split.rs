#![feature(test)]

extern crate plane_split;
extern crate test;

use std::sync::Arc;
use plane_split::{NaiveSplitter, Splitter, _make_grid};

#[bench]
fn bench_naive(b: &mut test::Bencher) {
    let polys = Arc::new(_make_grid(5));
    let mut splitter = NaiveSplitter::new();
    b.iter(|| {
        let p = polys.clone();
        splitter.solve(&p);
    });
}
