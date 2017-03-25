use std::collections::VecDeque;
use std::ops;
use {Polygon, Splitter};
use euclid::approxeq::ApproxEq;
use euclid::num::{One, Zero};

pub struct NaiveSplitter<T, U> {
    result: Vec<Polygon<T, U>>,
}

impl<
    T: Copy + PartialOrd + Zero + One + ApproxEq<T> +
       ops::Sub<T, Output=T> + ops::Add<T, Output=T> +
       ops::Mul<T, Output=T> + ops::Div<T, Output=T>,
    U,
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
                    let (res_add, res_add2) = existing.split(&line);
                    temp.push(res_add);
                    if let Some(res) = res_add2 {
                        temp.push(res);
                    }
                    let (new_todo, new_todo2) = polygon.split(&line);
                    heap_todo.push_back((new_todo, start_index + i + 1));
                    if let Some(todo) = new_todo2 {
                        heap_todo.push_back((todo, start_index + i + 1));
                    }
                }
            }
            self.result.extend(temp.drain(..));
        }
        &self.result
    }
}
