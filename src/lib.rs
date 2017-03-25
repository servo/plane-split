extern crate euclid;

mod naive;

use std::{mem, ops};
use euclid::{TypedPoint3D};
use euclid::approxeq::ApproxEq;
use euclid::num::{One, Zero};

pub use self::naive::NaiveSplitter;


pub type Index = u32;

pub struct Line<T, U> {
    pub origin: TypedPoint3D<T, U>,
    pub dir: TypedPoint3D<T, U>,
}

pub struct Polygon<T, U> {
    pub points: [TypedPoint3D<T, U>; 4],
    pub normal: TypedPoint3D<T, U>,
    pub offset: T,
    pub index: Index,
}

impl<T: Clone, U> Clone for Polygon<T, U> {
    fn clone(&self) -> Self {
        Polygon {
            points: [self.points[0].clone(),
                     self.points[1].clone(),
                     self.points[2].clone(),
                     self.points[3].clone()],
            normal: self.normal.clone(),
            offset: self.offset.clone(),
            index: self.index,
        }
    }
}

pub struct LineProjection<T> {
    pub markers: [T; 4],
}

impl<T: Copy + PartialOrd + ops::Sub<T, Output=T> + ops::Add<T, Output=T>> LineProjection<T> {
    pub fn get_bounds(&self) -> (T, T) {
        let (mut a, mut b, mut c, mut d) = (self.markers[0], self.markers[1], self.markers[2], self.markers[3]);
        // We could not just use `min/max` since they require `Ord` bound
        if a > c {
            mem::swap(&mut a, &mut c);
        }
        if b > d {
            mem::swap(&mut b, &mut d);
        }
        if a > b {
            mem::swap(&mut a, &mut b);
        }
        if c > d {
            mem::swap(&mut c, &mut d);
        }
        if b > c {
            mem::swap(&mut b, &mut c);
        }
        assert!(a <= b && b <= c && c <= d);
        (a, d)
    }

    pub fn intersect(&self, other: &Self) -> bool {
        // compute the bounds of both line projections
        let span = self.get_bounds();
        let other_span = other.get_bounds();
        // compute the total footprint
        let left = if span.0 < other_span.0 { span.0 } else { other_span.0 };
        let right = if span.1 > other_span.1 { span.1 } else { other_span.1 };
        // they intersect if the footprint is smaller than the sum
        right - left < span.1 - span.0 + other_span.1 - other_span.0
    }
}

fn scale<T: Copy + ops::Mul<T, Output=T>, U>(vec: TypedPoint3D<T, U>, factor: T) -> TypedPoint3D<T, U> {
    TypedPoint3D::new(vec.x * factor, vec.y * factor, vec.z * factor)
}

impl<T: Copy + PartialOrd + Zero + One + ApproxEq<T> +
        ops::Sub<T, Output=T> + ops::Add<T, Output=T> +
        ops::Mul<T, Output=T> + ops::Div<T, Output=T>,
     U> Polygon<T, U> {
    pub fn signed_distance_to(&self, point: &TypedPoint3D<T, U>) -> T {
        point.dot(self.normal) + self.offset
    }

    pub fn are_outside(&self, points: &[TypedPoint3D<T, U>]) -> bool {
        let d0 = self.signed_distance_to(&points[0]);
        points[1..].iter()
                   .find(|p| self.signed_distance_to(p) * d0 <= T::zero())
                   .is_none()
    }

    pub fn project_on(&self, vector: &TypedPoint3D<T, U>) -> LineProjection<T> {
        LineProjection {
            markers: [
                vector.dot(self.points[0]),
                vector.dot(self.points[1]),
                vector.dot(self.points[2]),
                vector.dot(self.points[3])
            ],
        }
    }

    pub fn intersect(&self, other: &Polygon<T, U>) -> Option<Line<T, U>> {
        if self.are_outside(&other.points) || other.are_outside(&self.points) {
            // one is completely outside the other
            return None
        }
        let cross_dir = self.normal.cross(other.normal);
        if cross_dir.dot(cross_dir) < T::approx_epsilon() {
            // polygons are co-planar
            return None
        }
        let self_proj = self.project_on(&cross_dir);
        let other_proj = other.project_on(&cross_dir);
        if !self_proj.intersect(&other_proj) {
            // projections on the line don't intersect
            return None
        }
        let center = scale(self.points[1] - self.points[0], other.offset) +
                     scale(other.points[1] - other.points[0], self.offset);
        Some(Line {
            origin: center,
            dir: cross_dir,
        })
    }

    pub fn split(&mut self, line: &Line<T, U>) -> (Polygon<T, U>, Option<Polygon<T, U>>) {
        let mut cuts = [None; 4];
        for ((&b, &a), cut) in self.points.iter()
                                        .cycle()
                                        .skip(1)
                                        .zip(self.points.iter())
                                        .zip(cuts.iter_mut()) {
            //a + (b-a) * t = line.origin + r * line.dir, where (t, r) = unknowns
            //a x line.dir + t * (b-a) x line.dir = line.origin x line.dir
            // t = (line.origin - a) x line.dir  / (b-a) x line.dir
            let c = b - a;
            let t = line.dir.cross(line.origin - a).x / line.dir.cross(c).x; //TODO
            if t > T::zero() && t < T::one() {
                *cut = Some(a + scale(c, t));
            }
        }
        let first = cuts.iter().position(|c| c.is_some()).unwrap();
        let second = first + cuts[first+1 ..].iter().position(|c| c.is_some()).unwrap();
        let (a, b) = (cuts[first].unwrap(), cuts[second].unwrap());
        match second-first {
            2 => {
                let mut other_points = self.points;
                other_points[first] = a;
                other_points[(first+3) % 4] = b;
                self.points[first+1] = a;
                self.points[first+2] = b;
                (Polygon {
                    points: other_points,
                    .. self.clone()
                }, None)
            }
            3 => {
                let xpoints = [
                    self.points[first+1],
                    self.points[(first+2) % 4],
                    self.points[(first+3) % 4],
                    b];
                let ypoints = [a, self.points[first+1], b, b];
                self.points = [self.points[first], a, b, b];
                (Polygon {
                    points: xpoints,
                    .. self.clone()
                }, Some(Polygon {
                    points: ypoints,
                    .. self.clone()
                }))
            }
            1 => {
                let xpoints = [
                    b,
                    self.points[(first+2) % 4],
                    self.points[(first+3) % 4],
                    self.points[first]
                    ];
                let ypoints = [self.points[first], a, b, b];
                self.points = [a, self.points[first+1], b, b];
                (Polygon {
                    points: xpoints,
                    .. self.clone()
                }, Some(Polygon {
                    points: ypoints,
                    .. self.clone()
                }))
            }
            _ => panic!("Unexpected indices {} {}", first, second),
        }
    }
}

pub trait Splitter<T, U> {
    fn solve(&mut self, polygons: &[Polygon<T, U>]) -> &[Polygon<T, U>];
}
