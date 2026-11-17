#include "meshcore/locate.hpp"

#include <cmath>
#include <stdexcept>

namespace meshcore {

Barycentric barycentric(Vec3 v0, Vec3 v1, Vec3 v2, Vec3 v3, Vec3 point) {
  const Vec3 a = v1 - v0;
  const Vec3 b = v2 - v0;
  const Vec3 c = v3 - v0;
  const Vec3 d = point - v0;
  Barycentric out;
  out.volume6 = scalar_triple(a, b, c);
  if (out.volume6 == 0.0) {
    return out;
  }
  out.w[1] = scalar_triple(d, b, c) / out.volume6;
  out.w[2] = scalar_triple(a, d, c) / out.volume6;
  out.w[3] = scalar_triple(a, b, d) / out.volume6;
  out.w[0] = 1.0 - out.w[1] - out.w[2] - out.w[3];
  return out;
}

bool is_inside(const Barycentric& b, double tol, bool require_positive_volume) {
  if (std::abs(b.volume6) <= 1.0e-18) {
    return false;
  }
  if (require_positive_volume && b.volume6 <= 0.0) {
    return false;
  }
  for (double weight : b.w) {
    if (weight < -tol) {
      return false;
    }
  }
  return true;
}

bool tet_contains(const Mesh& mesh, std::uint32_t tet, Vec3 point, double tol,
                  bool require_positive_volume) {
  const std::array<std::uint32_t, 4>& ids = mesh.tets[tet];
  return is_inside(barycentric(mesh.vertices[ids[0]], mesh.vertices[ids[1]], mesh.vertices[ids[2]],
                               mesh.vertices[ids[3]], point),
                   tol, require_positive_volume);
}

LocateHit locate_scan(const Mesh& mesh, Vec3 point, double tol, bool require_positive_volume) {
  LocateHit hit;
  const auto n = static_cast<std::uint32_t>(mesh.tets.size());
  for (std::uint32_t t = 0; t < n; ++t) {
    ++hit.barycentric_tests;
    if (tet_contains(mesh, t, point, tol, require_positive_volume)) {
      hit.tet = t;
      return hit;
    }
  }
  return hit;
}

LocateHit locate_walk(const Mesh& mesh, Vec3 point, std::uint32_t seed, double tol,
                      bool require_positive_volume) {
  if (mesh.neighbors.size() != mesh.tets.size()) {
    throw std::logic_error("locate_walk requires adjacency");
  }
  if (seed >= mesh.tets.size()) {
    throw std::invalid_argument("walk seed is out of range");
  }

  LocateHit hit;
  std::uint32_t current = seed;
  const std::uint32_t limit = static_cast<std::uint32_t>(mesh.tets.size()) + 1u;
  while (hit.steps <= limit) {
    ++hit.steps;
    const std::array<std::uint32_t, 4>& ids = mesh.tets[current];
    const Barycentric weights =
        barycentric(mesh.vertices[ids[0]], mesh.vertices[ids[1]], mesh.vertices[ids[2]],
                    mesh.vertices[ids[3]], point);
    ++hit.barycentric_tests;
    if (is_inside(weights, tol, require_positive_volume)) {
      hit.tet = current;
      return hit;
    }

    int worst = 0;
    for (int i = 1; i < 4; ++i) {
      if (weights.w[i] < weights.w[worst]) {
        worst = i;
      }
    }
    const std::int32_t next = mesh.neighbors[current][static_cast<std::size_t>(worst)];
    if (next < 0) {
      return hit;
    }
    current = static_cast<std::uint32_t>(next);
  }
  return hit;
}

}  // namespace meshcore
