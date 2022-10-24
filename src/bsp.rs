use crate::{Plane, PlaneCut, Polygon};

use euclid::default::{Point3D, Vector3D};

use std::fmt;

/// Binary Space Partitioning splitter, uses a BSP tree.
pub struct BspSplitter<A: Copy> {
    tree: BspNode<A>,
    result: Vec<Polygon<A>>,
}

impl<A: Copy> BspSplitter<A> {
    /// Create a new BSP splitter.
    pub fn new() -> Self {
        BspSplitter {
            tree: BspNode::new(),
            result: Vec::new(),
        }
    }
}

impl<A> BspSplitter<A>
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

    /// Process a set of polygons at once.
    pub fn solve(&mut self, input: &[Polygon<A>], view: Vector3D<f64>) -> &[Polygon<A>]
    where
        A: Copy,
    {
        self.reset();
        for p in input {
            self.add(p.clone());
        }
        self.sort(view)
    }
}

/// Add a list of planes to a particular front/end branch of some root node.
fn add_side<A: Copy>(side: &mut Option<Box<BspNode<A>>>, mut planes: Vec<Polygon<A>>) {
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
pub struct BspNode<A: Copy> {
    values: Vec<Polygon<A>>,
    front: Option<Box<BspNode<A>>>,
    back: Option<Box<BspNode<A>>>,
}

impl<A: Copy> BspNode<A> {
    /// Create a new node.
    pub fn new() -> Self {
        BspNode {
            values: Vec::new(),
            front: None,
            back: None,
        }
    }
}

impl<A: Copy> BspNode<A> {
    /// Insert a value into the sub-tree starting with this node.
    /// This operation may spawn additional leafs/branches of the tree.
    pub fn insert(&mut self, value: Polygon<A>) {
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
    pub fn order(&self, base: &Polygon<A>, out: &mut Vec<Polygon<A>>) {
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
