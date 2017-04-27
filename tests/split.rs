extern crate euclid;
extern crate plane_split;

use std::f32::consts::FRAC_PI_4;
use euclid::{Radians, TypedMatrix4D, TypedPoint2D, TypedPoint3D, TypedSize2D, TypedRect};
use plane_split::{BspSplitter, NaiveSplitter, Polygon, Splitter, _make_grid};


fn grid_impl(count: usize, splitter: &mut Splitter<f32, ()>) {
    let polys = _make_grid(count);
    let result = splitter.solve(&polys, TypedPoint3D::new(0.0, 0.0, 1.0));
    assert_eq!(result.len(), count + count*count + count*count*count);
}

#[test]
fn grid_naive() {
    grid_impl(2, &mut NaiveSplitter::new());
}

#[test]
fn grid_bsp() {
    grid_impl(2, &mut BspSplitter::new());
}


fn sort_impl(splitter: &mut Splitter<f32, ()>) {
    let transform0: TypedMatrix4D<f32, (), ()> =
        TypedMatrix4D::create_rotation(0.0, 1.0, 0.0, Radians::new(-FRAC_PI_4));
    let transform1: TypedMatrix4D<f32, (), ()> =
        TypedMatrix4D::create_rotation(0.0, 1.0, 0.0, Radians::new(0.0));
    let transform2: TypedMatrix4D<f32, (), ()> =
        TypedMatrix4D::create_rotation(0.0, 1.0, 0.0, Radians::new(FRAC_PI_4));

    let rect: TypedRect<f32, ()> = TypedRect::new(TypedPoint2D::new(-10.0, -10.0), TypedSize2D::new(20.0, 20.0));
    let polys = [
        Polygon::from_transformed_rect(rect, transform0, 0),
        Polygon::from_transformed_rect(rect, transform1, 1),
        Polygon::from_transformed_rect(rect, transform2, 2),
    ];

    let result = splitter.solve(&polys, TypedPoint3D::new(0.0, 0.0, -1.0));
    let ids: Vec<_> = result.iter().map(|poly| poly.anchor).collect();
    assert_eq!(&ids, &[2, 1, 0, 1, 2]);
}

#[test]
fn sort_naive() {
    sort_impl(&mut NaiveSplitter::new());
}

#[test]
fn sort_bsp() {
    sort_impl(&mut BspSplitter::new());
}
