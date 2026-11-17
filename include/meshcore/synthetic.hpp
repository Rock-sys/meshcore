#pragma once

#include "meshcore/vec3.hpp"

#include <cmath>

namespace meshcore {

// Regular tet, edge length 2√2, signed volume 8/3, circumradius √3.
// All three shape metrics below are exactly 1.
inline void regular_tet(Vec3& a, Vec3& b, Vec3& c, Vec3& d) {
  a = {1.0, 1.0, 1.0};
  b = {1.0, -1.0, -1.0};
  c = {-1.0, -1.0, 1.0};
  d = {-1.0, 1.0, -1.0};
}

// Equilateral base, apex 1e-4 above the centroid. Edge ratio stays under 2.
// Scaled Jacobian and radius ratio collapse. This is a sliver, not a needle.
inline void sliver_tet(Vec3& a, Vec3& b, Vec3& c, Vec3& d) {
  const double s = std::sqrt(3.0);
  a = {0.0, 0.0, 0.0};
  b = {1.0, 0.0, 0.0};
  c = {0.5, 0.5 * s, 0.0};
  d = {0.5, s / 6.0, 1.0e-4};
}

inline void inverted_tet(Vec3& a, Vec3& b, Vec3& c, Vec3& d) {
  regular_tet(a, b, c, d);
  const Vec3 swap = b;
  b = c;
  c = swap;
}

}  // namespace meshcore
