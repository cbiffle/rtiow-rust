use rand::prelude::*;
use std::ops::Range;

use crate::aabb::Aabb;
use crate::object::{HitRecord, Object};
use crate::ray::Ray;
use crate::vec3::Axis::{self, *};

#[derive(Debug)]
pub struct Bvh {
    bounding_box: Aabb,
    size: usize,
    contents: BvhContents,
}

#[derive(Debug)]
pub enum BvhContents {
    Node { left: Box<Bvh>, right: Box<Bvh> },
    Leaf(Box<dyn Object>),
}

impl Bvh {
    pub fn new(mut objs: Vec<Box<Object>>, exposure: Range<f32>, rng: &mut impl Rng) -> Self {
        // Note: though this BVH implementation is largely derived from Peter
        // Shirley's, it does *not* use the random axis selection and sort
        // routine, because it tended to fall into pathological cases.

        fn axis_range(objs: &[Box<Object>], exposure: Range<f32>, axis: Axis) -> f32 {
            let range = objs.iter().fold(std::f32::MAX..std::f32::MIN, |range, o| {
                let bb = o.bounding_box(exposure.clone());
                let min = bb.min[axis].min(bb.max[axis]);
                let max = bb.min[axis].max(bb.max[axis]);
                range.start.min(min)..range.end.max(max)
            });
            range.end - range.start
        }

        // Find the axis that has the greatest range for this set of objects.
        let axis = {
            let mut ranges = [
                (X, axis_range(&objs, exposure.clone(), X)),
                (Y, axis_range(&objs, exposure.clone(), Y)),
                (Z, axis_range(&objs, exposure.clone(), Z)),
            ];
            // Note reversed comparison function, to sort descending:
            ranges.sort_unstable_by(|a, b| b.1.partial_cmp(&a.1).unwrap());
            ranges[0].0
        };

        // Sort objects along it by centroid. (Actually, by centroid*2. This is
        // equivalent and cheaper.)
        objs.sort_unstable_by(|a, b| {
            let abb = a.bounding_box(exposure.clone());
            let bbb = b.bounding_box(exposure.clone());
            let av = abb.min[axis] + abb.max[axis];
            let bv = bbb.min[axis] + bbb.max[axis];
            av.partial_cmp(&bv).unwrap()
        });

        match objs.len() {
            0 => panic!("Can't create a BVH from zero objects."),
            1 => Bvh {
                bounding_box: objs[0].bounding_box(exposure.clone()),
                size: 1,
                contents: BvhContents::Leaf(objs.pop().unwrap()),
            },
            _ => {
                // Divide space at the median point of the selected axis.
                let right = Box::new(Bvh::new(
                    objs.drain(objs.len() / 2..).collect(),
                    exposure.clone(),
                    rng,
                ));
                let left = Box::new(Bvh::new(objs, exposure.clone(), rng));

                Bvh {
                    bounding_box: left.bounding_box.merge(right.bounding_box),
                    size: left.size + right.size,
                    contents: BvhContents::Node { left, right },
                }
            }
        }
    }
}

impl Object for Bvh {
    fn hit<'o>(&'o self, ray: &Ray, mut t_range: Range<f32>, rng: &mut dyn FnMut() -> f32) -> Option<HitRecord<'o>> {
        if self.bounding_box.hit(ray, t_range.clone()) {
            match &self.contents {
                BvhContents::Node { left, right } => {
                    let hit_left = left.hit(ray, t_range.clone(), rng);

                    // Don't bother searching past the left hit in the right
                    // space.
                    if let Some(h) = &hit_left {
                        t_range.end = h.t;
                    }

                    let hit_right = right.hit(ray, t_range, rng);

                    match (hit_left, hit_right) {
                        (h, None) | (None, h) => h,
                        (Some(hl), Some(hr)) => {
                            if hl.t < hr.t {
                                Some(hl)
                            } else {
                                Some(hr)
                            }
                        }
                    }
                }
                BvhContents::Leaf(obj) => obj.hit(ray, t_range, rng),
            }
        } else {
            None
        }
    }

    fn bounding_box(&self, exposure: std::ops::Range<f32>) -> Aabb {
        self.bounding_box
    }
}

// TODO this no longer has much value
pub fn from_scene(scene: Vec<Box<dyn Object>>, exposure: Range<f32>, rng: &mut impl Rng) -> Bvh {
    Bvh::new(scene, exposure, rng)
}
