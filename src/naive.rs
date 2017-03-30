use std::collections::VecDeque;
use std::{fmt, ops};
use {Polygon, Splitter};
use euclid::approxeq::ApproxEq;
use euclid::num::{One, Zero};

pub struct NaiveSplitter<T, U> {
    result: Vec<Polygon<T, U>>,
}

impl<T, U> NaiveSplitter<T, U> {
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

#[cfg(test)]
mod tests {
    use std::sync::Arc;
    use euclid::TypedPoint3D;
    #[cfg(test)]
    use test::Bencher;
    use super::*;

    fn make_parallel(count: usize) -> Vec<Polygon<f32, ()>> {
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

    #[test]
    fn naive() {
        let count = 2;
        let polys = make_parallel(count);
        let mut splitter = NaiveSplitter::new();
        let result = splitter.solve(&polys);
        assert_eq!(result.len(), 3*count*count*count);
    }

    #[bench]
    fn bench_naive(b: &mut Bencher) {
        let polys = Arc::new(make_parallel(5));
        let mut splitter = NaiveSplitter::new();
        b.iter(|| {
            let p = polys.clone();
            splitter.solve(&p);
        });
    }
}
