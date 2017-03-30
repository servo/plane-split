extern crate euclid;

mod naive;

use std::{fmt, mem, ops};
use euclid::TypedPoint3D;
use euclid::approxeq::ApproxEq;
use euclid::num::{One, Zero};

pub use self::naive::NaiveSplitter;

pub type Index = u32;

/// A generic line.
#[derive(Debug)]
pub struct Line<T, U> {
    /// Arbitrary point on the line.
    pub origin: TypedPoint3D<T, U>,
    /// Normalized direction of the line.
    pub dir: TypedPoint3D<T, U>,
}

impl<
    T: Copy + One + Zero + PartialEq + ApproxEq<T> + ops::Add<T, Output=T> + ops::Sub<T, Output=T> + ops::Mul<T, Output=T>,
    U,
> Line<T, U> {
    /// Check if the line has consistent parameters.
    pub fn is_valid(&self) -> bool {
        self.dir.dot(self.dir).approx_eq(&T::one())
    }
    /// Check if two lines match each other.
    pub fn matches(&self, other: &Self) -> bool {
        let diff = self.origin - other.origin;
        let zero = TypedPoint3D::new(T::zero(), T::zero(), T::zero());
        self.dir.cross(other.dir).approx_eq(&zero) &&
        self.dir.cross(diff).approx_eq(&zero)
    }
}

/// A convex flat polygon with 4 points, defined by equation:
/// dot(v, normal) + offset = 0
#[derive(Debug, PartialEq)]
pub struct Polygon<T, U> {
    /// Points making the polygon.
    pub points: [TypedPoint3D<T, U>; 4],
    /// Normalized vector perpendicular to the polygon plane.
    pub normal: TypedPoint3D<T, U>,
    /// Constant offset from the normal plane.
    pub offset: T,
    /// An original index of the polygon. Used to track the
    /// source when a polygon gets split.
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

/// The projection of a `Polygon` on a line.
pub struct LineProjection<T> {
    /// Projected value of each point in the polygon.
    pub markers: [T; 4],
}

impl<T: Copy + PartialOrd + ops::Sub<T, Output=T> + ops::Add<T, Output=T>> LineProjection<T> {
    pub fn get_bounds(&self) -> (T, T) {
        let (mut a, mut b, mut c, mut d) = (self.markers[0], self.markers[1], self.markers[2], self.markers[3]);
        // bitonic sort of 4 elements
        // we could not just use `min/max` since they require `Ord` bound
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
        debug_assert!(a <= b && b <= c && c <= d);
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

impl<T: Copy + fmt::Debug + PartialOrd + Zero + One + ApproxEq<T> +
        ops::Sub<T, Output=T> + ops::Add<T, Output=T> +
        ops::Mul<T, Output=T> + ops::Div<T, Output=T>,
     U> Polygon<T, U> {

    /// Return the signed distance from this polygon to a point.
    /// The distance is negative if the point is on the other side of the polygon
    /// from the direction of the normal.
    pub fn signed_distance_to(&self, point: &TypedPoint3D<T, U>) -> T {
        point.dot(self.normal) + self.offset
    }

    /// Check if all the points are indeed placed on the plane defined by
    /// the normal and offset, and the winding order is consistent.
    pub fn is_valid(&self) -> bool {
        let is_planar = self.points.iter()
                                   .find(|p| !self.signed_distance_to(p).approx_eq(&T::zero()))
                                   .is_none();
        let edges = [self.points[1] - self.points[0],
                     self.points[2] - self.points[1],
                     self.points[3] - self.points[2],
                     self.points[0] - self.points[3]];
        let anchor = edges[3].cross(edges[0]);
        let is_winding = edges.iter()
                              .zip(edges[1..].iter())
                              .find(|&(a, &b)| a.cross(b).dot(anchor) < T::zero())
                              .is_none();
        is_planar && is_winding
    }

    /// Check if a convex shape defined by a set of points is completely
    /// outside of this polygon. Merely touching the surface is not
    /// considered an intersection.
    pub fn are_outside(&self, points: &[TypedPoint3D<T, U>]) -> bool {
        let d0 = self.signed_distance_to(&points[0]);
        points[1..].iter()
                   .find(|p| self.signed_distance_to(p) * d0 <= T::zero())
                   .is_none()
    }

    /// Check if this polygon contains another one.
    pub fn contains(&self, other: &Self) -> bool {
        //TODO: actually check for inside/outside
        self.normal == other.normal && self.offset == other.offset
    }

    /// Project this polygon onto a 3D vector, returning a line projection.
    /// Note: we can think of it as a projection to a ray placed at the origin.
    pub fn project_on(&self, vector: &TypedPoint3D<T, U>) -> LineProjection<T> {
        LineProjection {
            markers: [
                vector.dot(self.points[0]),
                vector.dot(self.points[1]),
                vector.dot(self.points[2]),
                vector.dot(self.points[3]),
            ],
        }
    }

    pub fn intersect(&self, other: &Self) -> Option<Line<T, U>> {
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
        // compute any point on the intersection between planes
        // (n1, v) + d1 = 0
        // (n2, v) + d2 = 0
        // v = a*n1/w + b*n2/w; w = (n1, n2)
        // v = (d2*w - d1) / (1 - w*w) * n1 - (d2 - d1*w) / (1 - w*w) * n2
        let w = self.normal.dot(other.normal);
        let factor = T::one() / (T::one() - w * w);
        let center = scale(self.normal, (other.offset * w - self.offset) * factor) -
                     scale(other.normal, (other.offset - self.offset * w) * factor);
        Some(Line {
            origin: center,
            dir: cross_dir,
        })
    }

    pub fn split(&mut self, line: &Line<T, U>) -> (Option<Polygon<T, U>>, Option<Polygon<T, U>>) {
        // check if the cut is within the polygon plane first
        if !self.normal.dot(line.dir).approx_eq(&T::zero()) ||
           !self.signed_distance_to(&line.origin).approx_eq(&T::zero()) {
            return (None, None)
        }
        // compute the intersection points for each edge
        let mut cuts = [None; 4];
        for ((&b, &a), cut) in self.points.iter()
                                          .cycle()
                                          .skip(1)
                                          .zip(self.points.iter())
                                          .zip(cuts.iter_mut()) {
            // intersecting line segment [a, b] with `line`
            //a + (b-a) * t = r + k * d
            //(a, d) + t * (b-a, d) - (r, d) = k
            // a + t * (b-a) = r + t * (b-a, d) * d + (a-r, d) * d
            // t * ((b-a) - (b-a, d)*d) = (r-a) - (r-a, d) * d
            let pr = line.origin - a - scale(line.dir, line.dir.dot(line.origin - a));
            let pb = b - a - scale(line.dir, line.dir.dot(b - a));
            let denom = pb.dot(pb);
            if !denom.approx_eq(&T::zero()) {
                let t = pr.dot(pb) / denom;
                if t > T::zero() && t < T::one() {
                    *cut = Some(a + scale(b - a, t));
                }
            }
        }

        let first = match cuts.iter().position(|c| c.is_some()) {
            Some(pos) => pos,
            None => return (None, None),
        };
        let second = match cuts[first+1 ..].iter().position(|c| c.is_some()) {
            Some(pos) => first + 1 + pos,
            None => return (None, None),
        };
        //TODO: can be optimized for when the polygon has a redundant 4th vertex
        let (a, b) = (cuts[first].unwrap(), cuts[second].unwrap());
        match second-first {
            2 => {
                let mut other_points = self.points;
                other_points[first] = a;
                other_points[(first+3) % 4] = b;
                self.points[first+1] = a;
                self.points[first+2] = b;
                let poly = Polygon {
                    points: other_points,
                    .. self.clone()
                };
                (Some(poly), None)
            }
            3 => {
                let xpoints = [
                    self.points[first+1],
                    self.points[first+2],
                    self.points[first+3],
                    b];
                let ypoints = [a, self.points[first+1], b, b];
                self.points = [self.points[first], a, b, b];
                let poly1 = Polygon {
                    points: xpoints,
                    .. self.clone()
                };
                let poly2 = Polygon {
                    points: ypoints,
                    .. self.clone()
                };
                (Some(poly1), Some(poly2))
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
                let poly1 = Polygon {
                    points: xpoints,
                    .. self.clone()
                };
                let poly2 = Polygon {
                    points: ypoints,
                    .. self.clone()
                };
                (Some(poly1), Some(poly2))
            }
            _ => panic!("Unexpected indices {} {}", first, second),
        }
    }
}

pub trait Splitter<T, U> {
    fn solve(&mut self, polygons: &[Polygon<T, U>]) -> &[Polygon<T, U>];
}


#[doc(hidden)]
pub fn _make_grid(count: usize) -> Vec<Polygon<f32, ()>> {
    let mut polys: Vec<Polygon<f32, ()>> = Vec::with_capacity(count*3);
    let len = count as f32;
    polys.extend((0 .. count).map(|i| Polygon {
        points: [
            TypedPoint3D::new(0.0, i as f32, 0.0),
            TypedPoint3D::new(len, i as f32, 0.0),
            TypedPoint3D::new(len, i as f32, len),
            TypedPoint3D::new(0.0, i as f32, len),
        ],
        normal: TypedPoint3D::new(0.0, 1.0, 0.0),
        offset: -(i as f32),
        index: 1,
    }));
    polys.extend((0 .. count).map(|i| Polygon {
        points: [
            TypedPoint3D::new(i as f32, 0.0, 0.0),
            TypedPoint3D::new(i as f32, len, 0.0),
            TypedPoint3D::new(i as f32, len, len),
            TypedPoint3D::new(i as f32, 0.0, len),
        ],
        normal: TypedPoint3D::new(1.0, 0.0, 0.0),
        offset: -(i as f32),
        index: 1,
    }));
    polys.extend((0 .. count).map(|i| Polygon {
        points: [
            TypedPoint3D::new(0.0, 0.0, i as f32),
            TypedPoint3D::new(len, 0.0, i as f32),
            TypedPoint3D::new(len, len, i as f32),
            TypedPoint3D::new(0.0, len, i as f32),
        ],
        normal: TypedPoint3D::new(0.0, 0.0, 1.0),
        offset: -(i as f32),
        index: 1,
    }));
    polys
}
