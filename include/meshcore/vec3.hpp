#pragma once

#include <cmath>

namespace meshcore {

struct Vec3 {
  double x = 0;
  double y = 0;
  double z = 0;

  constexpr Vec3() = default;
  constexpr Vec3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}

  constexpr Vec3 operator+(Vec3 rhs) const { return {x + rhs.x, y + rhs.y, z + rhs.z}; }
  constexpr Vec3 operator-(Vec3 rhs) const { return {x - rhs.x, y - rhs.y, z - rhs.z}; }
  constexpr Vec3 operator*(double s) const { return {x * s, y * s, z * s}; }
  constexpr Vec3 operator/(double s) const { return {x / s, y / s, z / s}; }
};

constexpr Vec3 operator*(double s, Vec3 v) { return v * s; }

constexpr double dot(Vec3 a, Vec3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }

constexpr Vec3 cross(Vec3 a, Vec3 b) {
  return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}

inline double length(Vec3 a) { return std::sqrt(dot(a, a)); }

constexpr double scalar_triple(Vec3 a, Vec3 b, Vec3 c) { return dot(a, cross(b, c)); }

inline double coord(Vec3 v, int axis) {
  if (axis == 0) {
    return v.x;
  }
  if (axis == 1) {
    return v.y;
  }
  return v.z;
}

}  // namespace meshcore
