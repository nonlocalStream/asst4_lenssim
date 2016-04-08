#include "bvh.h"

#include "CGL/CGL.h"
#include "static_scene/triangle.h"

#include <iostream>
#include <stack>

using namespace std;

namespace CGL { namespace StaticScene {

BVHAccel::BVHAccel(const std::vector<Primitive *> &_primitives,
                   size_t max_leaf_size) {

  root = construct_bvh(_primitives, max_leaf_size);

}

BVHAccel::~BVHAccel() {
  if (root) delete root;
}

BBox BVHAccel::get_bbox() const {
  return root->bb;
}

void BVHAccel::draw(BVHNode *node, const Color& c) const {
  if (node->isLeaf()) {
    for (Primitive *p : *(node->prims))
      p->draw(c);
  } else {
    draw(node->l, c);
    draw(node->r, c);
  }
}

void BVHAccel::drawOutline(BVHNode *node, const Color& c) const {
  if (node->isLeaf()) {
    for (Primitive *p : *(node->prims))
      p->drawOutline(c);
  } else {
    drawOutline(node->l, c);
    drawOutline(node->r, c);
  }
}

int BVHAccel::largest_axis(const Vector3D extent) {
  if (extent[0]>=extent[2]) {
    if (extent[0]>=extent[1])
        return 0;
    else
        return 1;
  }
  return 2;
}
BVHNode *BVHAccel::construct_bvh(const std::vector<Primitive*>& prims, size_t max_leaf_size) {
  
  // TODO Part 2, task 1:
  // Construct a BVH from the given vector of primitives and maximum leaf
  // size configuration. The starter code build a BVH aggregate with a
  // single leaf node (which is also the root) that encloses all the
  // primitives.

  BBox centroid_box, bbox;

  for (Primitive *p : prims) {
    BBox bb = p->get_bbox();
    bbox.expand(bb);
    Vector3D c = bb.centroid();
    centroid_box.expand(c);
  }

  // You'll want to adjust this code.
  // Right now we just return a single node containing all primitives.
  BVHNode *node = new BVHNode(bbox);
  node->bb = bbox;
  node->l = NULL;
  node->r = NULL;
  if (prims.size() <= max_leaf_size) {
    node->prims = new vector<Primitive *>(prims);
  } else {
    int axis = largest_axis(bbox.extent);
    double split = ((centroid_box.max+centroid_box.min)/2)[axis];
    std::vector<Primitive*> lprims, rprims, rest;
    for (Primitive *p : prims) {
      BBox bb = p->get_bbox();
      Vector3D c = bb.centroid();
      if (c[axis] < split) {
        lprims.push_back(p);
      } else if (c[axis] > split){
        rprims.push_back(p);
      } else {
        rest.push_back(p);
      }
    }
    for (int i = 0; i < rest.size(); i++) {
      if (lprims.size() < rprims.size()) {
        lprims.push_back(rest[i]);
      } else {
        rprims.push_back(rest[i]);  
      }
    }
    if (!lprims.empty()) {
      BVHNode* leftn = construct_bvh(lprims, max_leaf_size);
      node->l = leftn;
    }
    if (!rprims.empty()) {
      BVHNode* rightn = construct_bvh(rprims, max_leaf_size);
      node->r = rightn;
    }
  }
  return node;
}


bool BVHAccel::intersect(const Ray& ray, BVHNode *node) const {

  // TODO Part 2, task 3:
  // Implement BVH intersection.
  // Currently, we just naively loop over every primitive.
  if (node->isLeaf()) {
    for (Primitive *p : *(node->prims)) {
      total_isects++;
      if (p->intersect(ray)) 
        return true;
    }
  } else {
    double tmin = ray.min_t;
    double tmax = ray.max_t;
    if ((node->l != NULL) && (node->l->bb).intersect(ray, tmin, tmax)) {
        if (intersect(ray, node->l)) return true;
    };
    tmin = ray.min_t;
    tmax = ray.max_t;
    if ((node->r != NULL) && (node->r->bb).intersect(ray, tmin, tmax)) {
        if (intersect(ray, node->r)) return true;
    };
  }
  /**
  double tmin = ray.min_t;
  double tmax = ray.max_t;
  if ((root->bb).intersect(ray, tmin, tmax)) {
      cout << "intersected root" << endl;
  }
  */
  return false;
}

bool BVHAccel::intersect(const Ray& ray, Intersection* i, BVHNode *node) const {

  // TODO Part 2, task 3:
  // Implement BVH intersection.
  // Currently, we just naively loop over every primitive.
  bool hit = false;
  if (node->isLeaf()) {
    for (Primitive *p : *(node->prims)) {
      total_isects++;
      if (p->intersect(ray, i)) 
        hit = true;
    }
  } else {
    double tmin = ray.min_t;
    double tmax = ray.max_t;
    if ((node->l != NULL) && (node->l->bb).intersect(ray, tmin, tmax)) {
      if (intersect(ray, i, node->l)) {
        hit = true;
      };
    };
    tmin = ray.min_t;
    tmax = ray.max_t;
    if ((node->r != NULL) && (node->r->bb).intersect(ray, tmin, tmax)) {
      if (intersect(ray, i, node->r)) {
        hit = true;
      }
    };
  }
  /**
  bool hit = false;
  double tmin = ray.min_t;
  double tmax = ray.max_t;
  if ((root->bb).intersect(ray, tmin, tmax)) {
      cout << "intersected root" << endl;
  }
*/
  return hit;
}

}  // namespace StaticScene
}  // namespace CGL
