use std::{fmt, ops};
use {Polygon, Splitter};
use euclid::TypedPoint3D;
use euclid::approxeq::ApproxEq;
use num_traits::{One, Zero};


/// Naive plane splitter, has at least O(n^2) complexity.
pub struct NaiveSplitter<T, U> {
    result: Vec<Polygon<T, U>>,
    current: Vec<Polygon<T, U>>,
    temp: Vec<Polygon<T, U>>,
}

impl<T, U> NaiveSplitter<T, U> {
    /// Create a new `NaiveSplitter`.
    pub fn new() -> Self {
        NaiveSplitter {
            result: Vec::new(),
            current: Vec::new(),
            temp: Vec::new(),
        }
    }
}

impl<
    T: Copy + fmt::Debug + PartialOrd + Zero + One + ApproxEq<T> +
       ops::Sub<T, Output=T> + ops::Add<T, Output=T> +
       ops::Mul<T, Output=T> + ops::Div<T, Output=T>,
    U: fmt::Debug,
> Splitter<T, U> for NaiveSplitter<T, U> {
    fn reset(&mut self) {
        self.result.clear();
        self.current.clear();
        self.temp.clear();
    }

    fn get_all(&self) -> &[Polygon<T ,U>] {
        &self.result
    }

    fn add(&mut self, poly: Polygon<T, U>) -> &[Polygon<T, U>] {
        // "current" accumulates all the subdivisions of the originally
        // added polygon
        self.current.push(poly);
        for old in self.result.iter() {
            for new in self.current.iter_mut() {
                // temp accumulates all the new subdivisions to be added
                // to the current, since we can't modify it in place
                if let Some(line) = old.intersect(new) {
                    let (res_add1, res_add2) = new.split(&line);
                    if let Some(res) = res_add1 {
                        self.temp.push(res);
                    }
                    if let Some(res) = res_add2 {
                        self.temp.push(res);
                    }
                }
            }
            self.current.extend(self.temp.drain(..));
        }
        let index = self.result.len();
        self.result.extend(self.current.drain(..));
        &self.result[index..]
    }

    fn sort(&mut self, view: TypedPoint3D<T, U>) {
        //unimplemented!()
    }
}
