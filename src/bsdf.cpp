#include "bsdf.h"

#include <iostream>
#include <algorithm>
#include <utility>

using std::min;
using std::max;
using std::swap;

namespace CGL {

void make_coord_space(Matrix3x3& o2w, const Vector3D& n) {

    Vector3D z = Vector3D(n.x, n.y, n.z);
    Vector3D h = z;
    if (fabs(h.x) <= fabs(h.y) && fabs(h.x) <= fabs(h.z)) h.x = 1.0;
    else if (fabs(h.y) <= fabs(h.x) && fabs(h.y) <= fabs(h.z)) h.y = 1.0;
    else h.z = 1.0;

    z.normalize();
    Vector3D y = cross(h, z);
    y.normalize();
    Vector3D x = cross(z, y);
    x.normalize();

    o2w[0] = x;
    o2w[1] = y;
    o2w[2] = z;
}

// Diffuse BSDF //

Spectrum DiffuseBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  return albedo * (1.0 / PI);
}

Spectrum DiffuseBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
  *wi = sampler.get_sample(pdf);
  return albedo * (1.0 / PI);
}

// Mirror BSDF //

Spectrum MirrorBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  Vector3D refl;
  reflect(wo, &refl);
  double cos_w = dot(Vector3D(0,0,1), wo.unit());
  if (refl == wi) {
    return reflectance / cos_w;
  } else {
    return Spectrum();
  }
}

Spectrum MirrorBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {

  // TODO Part 5:
  // Implement MirrorBSDF
  Vector3D refl_in;
  reflect(wo, wi);
  double cos_w = dot(Vector3D(0,0,1), wo.unit());
  *pdf = 1.f;
  return reflectance / cos_w;
}

// Glossy BSDF //

/*
Spectrum GlossyBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  return Spectrum();
}

Spectrum GlossyBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
  *pdf = 1.0f;
  return reflect(wo, wi, reflectance);
}
*/

// Refraction BSDF //

Spectrum RefractionBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  return Spectrum();
}

Spectrum RefractionBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {

  // TODO Part 5:
  // Implement RefractionBSDF
  return Spectrum();
}

// Glass BSDF //

Spectrum GlassBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  return Spectrum();
}

Spectrum GlassBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {

  // TODO Part 5:
  // Compute Fresnel coefficient and either reflect or refract based on it.
  Vector3D refl, refr;
  reflect(wo, &refl);
  bool has_refr = refract(wo, &refr, ior);
  double cos_o = abs_cos_theta(wo.unit());
  if (!has_refr) {
      *wi = refl;
      *pdf = 1.f;
      return reflectance / cos_o;
  } else {
      double R0 = pow((1-ior)/(1+ior),2);
      double cos_theta = abs_cos_theta(refl.unit());
      double R = clamp(R0 + (1-R0)*pow((1-cos_theta),5), 0, 1);
      if (coin_flip(R)) {
        *wi = refl;
        *pdf = R;
        return R * reflectance / cos_o;
      } else {
        *wi = refr;
        *pdf = 1.f - R;
        double ni_div_no = (wo.z > 0) ? ior: 1/ior;
        return (1.f-R)*transmittance*pow(ni_div_no, 2)/cos_o;
      }
  }
  return Spectrum();
}

void BSDF::reflect(const Vector3D& wo, Vector3D* wi) {

  // TODO Part 5:
  // Implement reflection of wo about normal (0,0,1) and store result in wi.
  *wi = Vector3D(-wo.x, -wo.y, wo.z);
}

bool BSDF::refract(const Vector3D& wo, Vector3D* wi, float ior) {

  // TODO Part 5:
  // Use Snell's Law to refract wo surface and store result ray in wi.
  // Return false if refraction does not occur due to total internal reflection
  // and true otherwise. When dot(wo,n) is positive, then wo corresponds to a
  // ray entering the surface through vacuum.
  Vector3D n = Vector3D(0,0,1); 
  float ni, no, sin_theta_i, sin_theta_o;
  if (dot(wo,n) > 0) { //what about == 0?
    no = 1;
    ni = ior;
  } else {
    no = ior;
    ni = 1;
  }
  sin_theta_o = sin_theta(wo);
  if (no*sin_theta_o >= ni) return false;
  sin_theta_i = no*sin_theta_o/ni;
  double x2_plus_y2 = pow(wo.x,2)+pow(wo.y,2);
  double z2 = x2_plus_y2/pow(sin_theta_i,2) - x2_plus_y2;
  double z = sqrt(z2);

  if (dot(wo,n) > 0) { //if wo from above, z should be negative
      z = -z;
  }
  *wi = Vector3D(-wo.x, -wo.y, z).unit();
  return true;
}

// Emission BSDF //

Spectrum EmissionBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  return Spectrum();
}

Spectrum EmissionBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
  *pdf = 1.0 / PI;
  *wi  = sampler.get_sample(pdf);
  return Spectrum();
}

} // namespace CGL
