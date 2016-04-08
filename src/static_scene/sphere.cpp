#include "sphere.h"

#include <cmath>
#include <cfloat>

#include  "../bsdf.h"
#include "../misc/sphere_drawing.h"

namespace CGL { namespace StaticScene {

bool Sphere::test(const Ray& r, double& t1, double& t2) const {

  // TODO Part 1, task 4:
  // Implement ray - sphere intersection test.
  // Return true if there are intersections and writing the
  // smaller of the two intersection times in t1 and the larger in t2.
  double a = dot(r.d, r.d);
  double b = 2*dot(r.o-o, r.d);
  double c = dot(r.o-o, r.o-o) - pow(this->r,2);
  if (pow(b,2) - 4*a*c < -DBL_EPSILON) return false;
  if (a == 0) return false;
  double delta = pow(b,2) - 4*a*c;
  if (delta < 0) {
    delta = 0;
  }
  t1 = (-b-sqrt(delta))/2/a;
  t2 = (-b+sqrt(delta))/2/a;
  return true;
}

bool Sphere::intersect(const Ray& r) const {

  // TODO Part 1, task 4:
  // Implement ray - sphere intersection.
  // Note that you might want to use the the Sphere::test helper here.
  double t1, t2;
  if (test(r, t1, t2)) {
      if ((t1>=0) && (t1 >= r.min_t) && (t1 <= r.max_t)) {
        r.max_t = t1;
      } else if ((t2>=0)&&(t2 >= r.min_t) && (t2 <= r.max_t)) {
        r.max_t = t2;
      } else {
        return false;
      }
      return true;
  }
  return false;
}

bool Sphere::intersect(const Ray& r, Intersection *i) const {

  // TODO Part 1m task 4:
  // Implement ray - sphere intersection.
  // Note again that you might want to use the the Sphere::test helper here.
  // When an intersection takes place, the Intersection data should be updated
  // correspondingly.
  double t1, t2;
  if (test(r, t1, t2)) {
      if ((t1>=0) && (t1 >= r.min_t) && (t1 <= r.max_t)) {
        r.max_t = t1;
        i->t = t1;
        i->n = (r.o+t1*r.d-o).unit();
      } else if ((t2>=0) && (t2 >= r.min_t) && (t2 <= r.max_t)) {
        r.max_t = t2;
        i->t = t2;
        i->n = (r.o+t2*r.d-o).unit();
      } else {
        return false;
      }
      i->primitive = this;
      i->bsdf = get_bsdf();
      return true;
  }
  return false;
}

void Sphere::draw(const Color& c) const {
  Misc::draw_sphere_opengl(o, r, c);
}

void Sphere::drawOutline(const Color& c) const {
    //Misc::draw_sphere_opengl(o, r, c);
}


} // namespace StaticScene
} // namespace CGL
