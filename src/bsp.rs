use crate::{is_zero, Intersection, Plane, Polygon, Splitter};

use euclid::default::{Point3D, Vector3D};

use std::{fmt, iter};

/// A plane abstracted to the matter of partitioning.
pub trait BspPlane: Sized + Clone {
    /// Try to cut a different plane by this one.
    fn cut(&self, other: Self) -> PlaneCut<Self>;
    /// Check if a different plane is aligned in the same direction
    /// as this one.
    fn is_aligned(&self, otgher: &Self) -> bool;
}

impl<A> BspPlane for Polygon<A>
where
    A: Copy + fmt::Debug,
{
    fn cut(&self, mut poly: Self) -> PlaneCut<Self> {
        log::debug!("\tCutting anchor {:?} by {:?}", poly.anchor, self.anchor);
        log::trace!("\t\tbase {:?}", self.plane);

        //Note: we treat `self` as a plane, and `poly` as a concrete polygon here
        let (intersection, dist) = match self.plane.intersect(&poly.plane) {
            None => {
                let ndot = self.plane.normal.dot(poly.plane.normal);
                log::debug!("\t\tNormals are aligned with {:?}", ndot);
                let dist = self.plane.offset - ndot * poly.plane.offset;
                (Intersection::Coplanar, dist)
            }
            Some(_) if self.plane.are_outside(&poly.points[..]) => {
                //Note: we can't start with `are_outside` because it's subject to FP precision
                let dist = self.plane.signed_distance_sum_to(&poly);
                (Intersection::Outside, dist)
            }
            Some(line) => {
                //Note: distance isn't relevant here
                (Intersection::Inside(line), 0.0)
            }
        };

        match intersection {
            //Note: we deliberately make the comparison wider than just with T::epsilon().
            // This is done to avoid mistakenly ordering items that should be on the same
            // plane but end up slightly different due to the floating point precision.
            Intersection::Coplanar if is_zero(dist) => {
                log::debug!("\t\tCoplanar at {:?}", dist);
                PlaneCut::Sibling(poly)
            }
            Intersection::Coplanar | Intersection::Outside => {
                log::debug!("\t\tOutside at {:?}", dist);
                if dist > 0.0 {
                    PlaneCut::Cut {
                        front: vec![poly],
                        back: vec![],
                    }
                } else {
                    PlaneCut::Cut {
                        front: vec![],
                        back: vec![poly],
                    }
                }
            }
            Intersection::Inside(line) => {
                log::debug!("\t\tCut across {:?}", line);
                let (res_add1, res_add2) = poly.split_with_normal(&line, &self.plane.normal);
                let mut front = Vec::new();
                let mut back = Vec::new();

                for sub in iter::once(poly)
                    .chain(res_add1)
                    .chain(res_add2)
                    .filter(|p| !p.is_empty())
                {
                    let dist = self.plane.signed_distance_sum_to(&sub);
                    if dist > 0.0 {
                        log::trace!("\t\t\tdist {:?} -> front: {:?}", dist, sub);
                        front.push(sub)
                    } else {
                        log::trace!("\t\t\tdist {:?} -> back: {:?}", dist, sub);
                        back.push(sub)
                    }
                }

                PlaneCut::Cut { front, back }
            }
        }
    }

    fn is_aligned(&self, other: &Self) -> bool {
        self.plane.normal.dot(other.plane.normal) > 0.0
    }
}

/// Binary Space Partitioning splitter, uses a BSP tree.
pub struct BspSplitter<A> {
    tree: BspNode<Polygon<A>>,
    result: Vec<Polygon<A>>,
}

impl<A> BspSplitter<A> {
    /// Create a new BSP splitter.
    pub fn new() -> Self {
        BspSplitter {
            tree: BspNode::new(),
            result: Vec::new(),
        }
    }
}

impl<A> Splitter<A> for BspSplitter<A>
where
    A: Copy + fmt::Debug + Default,
{
    fn reset(&mut self) {
        self.tree = BspNode::new();
    }

    fn add(&mut self, poly: Polygon<A>) {
        self.tree.insert(poly);
    }

    fn sort(&mut self, view: Vector3D<f64>) -> &[Polygon<A>] {
        //debug!("\t\ttree before sorting {:?}", self.tree);
        let poly = Polygon {
            points: [Point3D::origin(); 4],
            plane: Plane {
                normal: -view, //Note: BSP `order()` is back to front
                offset: 0.0,
            },
            anchor: A::default(),
        };
        self.result.clear();
        self.tree.order(&poly, &mut self.result);
        &self.result
    }
}

/// The result of one plane being cut by another one.
/// The "cut" here is an attempt to classify a plane as being
/// in front or in the back of another one.
#[derive(Debug)]
pub enum PlaneCut<T> {
    /// The planes are one the same geometrical plane.
    Sibling(T),
    /// Planes are different, thus we can either determine that
    /// our plane is completely in front/back of another one,
    /// or split it into these sub-groups.
    Cut {
        /// Sub-planes in front of the base plane.
        front: Vec<T>,
        /// Sub-planes in the back of the base plane.
        back: Vec<T>,
    },
}

/// Add a list of planes to a particular front/end branch of some root node.
fn add_side<T: BspPlane>(side: &mut Option<Box<BspNode<T>>>, mut planes: Vec<T>) {
    if planes.len() != 0 {
        if side.is_none() {
            *side = Some(Box::new(BspNode::new()));
        }
        let node = side.as_mut().unwrap();
        for p in planes.drain(..) {
            node.insert(p)
        }
    }
}

/// A node in the `BspTree`, which can be considered a tree itself.
#[derive(Clone, Debug)]
pub struct BspNode<T> {
    values: Vec<T>,
    front: Option<Box<BspNode<T>>>,
    back: Option<Box<BspNode<T>>>,
}

impl<T> BspNode<T> {
    /// Create a new node.
    pub fn new() -> Self {
        BspNode {
            values: Vec::new(),
            front: None,
            back: None,
        }
    }
}

impl<T: BspPlane> BspNode<T> {
    /// Insert a value into the sub-tree starting with this node.
    /// This operation may spawn additional leafs/branches of the tree.
    pub fn insert(&mut self, value: T) {
        if self.values.is_empty() {
            self.values.push(value);
            return;
        }
        match self.values[0].cut(value) {
            PlaneCut::Sibling(value) => self.values.push(value),
            PlaneCut::Cut { front, back } => {
                add_side(&mut self.front, front);
                add_side(&mut self.back, back);
            }
        }
    }

    /// Build the draw order of this sub-tree into an `out` vector,
    /// so that the contained planes are sorted back to front according
    /// to the view vector defined as the `base` plane front direction.
    pub fn order(&self, base: &T, out: &mut Vec<T>) {
        let (former, latter) = match self.values.first() {
            None => return,
            Some(ref first) if base.is_aligned(first) => (self.front.as_ref(), self.back.as_ref()),
            Some(_) => (self.back.as_ref(), self.front.as_ref()),
        };

        if let Some(ref node) = former {
            node.order(base, out);
        }

        out.extend_from_slice(&self.values);

        if let Some(ref node) = latter {
            node.order(base, out);
        }
    }
}
