use std::collections::VecDeque;
use std::{fmt, ops};
use {Polygon, Splitter};
use euclid::approxeq::ApproxEq;
use euclid::num::{One, Zero};


/// Naive plane splitter, has at least O(n^2) complexity.
pub struct NaiveSplitter<T, U> {
    result: Vec<Polygon<T, U>>,
}

impl<T, U> NaiveSplitter<T, U> {
    /// Create a new `NaiveSplitter`.
    pub fn new() -> Self {
        NaiveSplitter {
            result: Vec::new(),
        }
    }
}

impl<
    T: Copy + fmt::Debug + PartialOrd + Zero + One + ApproxEq<T> +
       ops::Sub<T, Output=T> + ops::Add<T, Output=T> +
       ops::Mul<T, Output=T> + ops::Div<T, Output=T>,
    U: fmt::Debug,
> Splitter<T, U> for NaiveSplitter<T, U> {
    fn solve(&mut self, polygons: &[Polygon<T, U>]) -> &[Polygon<T, U>] {
        self.result.clear();
        let mut temp = Vec::new();
        let mut heap_todo: VecDeque<_> = polygons.iter()
                                                 .map(|p| (p.clone(), 0usize))
                                                 .collect();
        while let Some((mut polygon, start_index)) = heap_todo.pop_front() {
            for (i, existing) in self.result[start_index..].iter_mut().enumerate() {
                if let Some(line) = polygon.intersect(existing) {
                    let (res_add1, res_add2) = existing.split(&line);
                    if let Some(res) = res_add1 {
                        temp.push(res);
                    }
                    if let Some(res) = res_add2 {
                        temp.push(res);
                    }
                    let (new_todo1, new_todo2) = polygon.split(&line);
                    if let Some(todo) = new_todo1 {
                        heap_todo.push_back((todo, start_index + i + 1));
                    }
                    if let Some(todo) = new_todo2 {
                        heap_todo.push_back((todo, start_index + i + 1));
                    }
                }
            }
            self.result.push(polygon);
            self.result.extend(temp.drain(..));
        }
        &self.result
    }
}
