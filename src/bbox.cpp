#include "bbox.h"

#include "GL/glew.h"

#include <algorithm>
#include <iostream>

namespace CGL {

void BBox::intersect_planes(const Ray& r, int axis, double& tmin, double& tmax) const{
    double p1 = min[axis];
    double p2 = max[axis];
    double o = r.o[axis];
    double d = r.d[axis];
    if (d == 0.f) {
      if ((p1<=o) && (o<=p2)) {
        tmin = -INF_D;
        tmax = INF_D;
      } else {
        tmin = INF_D;
        tmax = -INF_D;
      }
    } else {
      double t1 = (p1-o)/d;
      double t2 = (p2-o)/d;
      tmin = std::min(t1,t2);
      tmax = std::max(t1,t2);
    }
}
bool BBox::intersect(const Ray& r, double& t0, double& t1) const {
  // TODO Part 2, task 2:
  // Implement ray - bounding box intersection test
  // If the ray intersected the bounding box within the range given by
  // t0, t1, update t0 and t1 with the new intersection times.
  //std::cout << "testing bbox: " << std::endl;
  double tminx, tminy, tminz;
  double tmaxx, tmaxy, tmaxz;
  intersect_planes(r, 0, tminx, tmaxx);
  intersect_planes(r, 1, tminy, tmaxy);
  intersect_planes(r, 2, tminz, tmaxz);
  double tmin = std::max(std::max(tminx, tminy), tminz);
  double tmax = std::min(std::min(tmaxx, tmaxy), tmaxz);
  if ((tmin <= tmax)) {//&&(t0<=tmin)&&(t1>=tmax)) {
      t0 = tmin;
      t1 = tmax;
      return true;
  }
  return false;
/*
    Vector3D invdir = 1 / r.d; 
    sign[0] = (invdir.x < 0); 
    sign[1] = (invdir.y < 0); 
    sign[2] = (invdir.z < 0); 
    float tmin, tmax, tymin, tymax, tzmin, tzmax; 
 
    tmin = (bounds[r.sign[0]].x - r.orig.x) * r.invdir.x; 
    tmax = (bounds[1-r.sign[0]].x - r.orig.x) * r.invdir.x; 
    tymin = (bounds[r.sign[1]].y - r.orig.y) * r.invdir.y; 
    tymax = (bounds[1-r.sign[1]].y - r.orig.y) * r.invdir.y; 
 
    if ((tmin > tymax) || (tymin > tmax)) 
        return false; 
    if (tymin > tmin) 
        tmin = tymin; 
    if (tymax < tmax) 
        tmax = tymax; 
 
    tzmin = (bounds[r.sign[2]].z - r.orig.z) * r.invdir.z; 
    tzmax = (bounds[1-r.sign[2]].z - r.orig.z) * r.invdir.z; 
 
    if ((tmin > tzmax) || (tzmin > tmax)) 
        return false; 
    if (tzmin > tmin) 
        tmin = tzmin; 
    if (tzmax < tmax) 
        tmax = tzmax; 
  */
 
    return true;     
}

void BBox::draw(Color c) const {

  glColor4f(c.r, c.g, c.b, c.a);

	// top
	glBegin(GL_LINE_STRIP);
	glVertex3d(max.x, max.y, max.z);
  glVertex3d(max.x, max.y, min.z);
  glVertex3d(min.x, max.y, min.z);
  glVertex3d(min.x, max.y, max.z);
  glVertex3d(max.x, max.y, max.z);
	glEnd();

	// bottom
	glBegin(GL_LINE_STRIP);
  glVertex3d(min.x, min.y, min.z);
  glVertex3d(min.x, min.y, max.z);
  glVertex3d(max.x, min.y, max.z);
  glVertex3d(max.x, min.y, min.z);
  glVertex3d(min.x, min.y, min.z);
	glEnd();

	// side
	glBegin(GL_LINES);
	glVertex3d(max.x, max.y, max.z);
  glVertex3d(max.x, min.y, max.z);
	glVertex3d(max.x, max.y, min.z);
  glVertex3d(max.x, min.y, min.z);
	glVertex3d(min.x, max.y, min.z);
  glVertex3d(min.x, min.y, min.z);
	glVertex3d(min.x, max.y, max.z);
  glVertex3d(min.x, min.y, max.z);
	glEnd();

}

std::ostream& operator<<(std::ostream& os, const BBox& b) {
  return os << "BBOX(" << b.min << ", " << b.max << ")";
}

} // namespace CGL
