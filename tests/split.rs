extern crate euclid;
extern crate plane_split;

use std::f32::consts::FRAC_PI_4;
use euclid::{Radians, TypedMatrix4D, TypedPoint2D, TypedPoint3D, TypedSize2D, TypedRect};
use plane_split::{NaiveSplitter, Polygon, Splitter, _make_grid};

#[test]
fn naive_grid() {
    let count = 2;
    let polys = _make_grid(count);
    let mut splitter = NaiveSplitter::new();
    for poly in polys {
        splitter.add(poly);
    }
    assert_eq!(splitter.get_all().len(), count + count*count + count*count*count);
}

#[test]
fn naive_sort() {
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

    let mut splitter = NaiveSplitter::new();
    splitter.solve(&polys, TypedPoint3D::new(0.0, 0.0, -1.0));
    let ids: Vec<_> = splitter.get_all().iter().map(|poly| poly.anchor).collect();
    assert_eq!(&ids, &[2, 1, 0, 1, 2]);
}
