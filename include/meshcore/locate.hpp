#pragma once

#include "meshcore/mesh.hpp"
#include "meshcore/vec3.hpp"

#include <cstdint>
#include <optional>

namespace meshcore {

struct Barycentric {
  double w[4] = {0, 0, 0, 0};
  double volume6 = 0;
};

struct LocateHit {
  std::optional<std::uint32_t> tet;
  std::uint32_t barycentric_tests = 0;
  std::uint32_t nodes_visited = 0;
  std::uint32_t steps = 0;
};

[[nodiscard]] Barycentric barycentric(Vec3 v0, Vec3 v1, Vec3 v2, Vec3 v3, Vec3 point);

[[nodiscard]] bool is_inside(const Barycentric& b, double tol, bool require_positive_volume);

[[nodiscard]] bool tet_contains(const Mesh& mesh, std::uint32_t tet, Vec3 point, double tol,
                                bool require_positive_volume);

// require_positive_volume rejects inverted elements. Location and orientation
// are separate checks: a negative Jacobian still has a geometric interior.
[[nodiscard]] LocateHit locate_scan(const Mesh& mesh, Vec3 point, double tol,
                                    bool require_positive_volume);

// Straight walk across the face of the most negative barycentric coordinate.
// The seed must lie in the same connected component. A hop limit stops cycles
// when a sliver makes the orientation test unreliable.
[[nodiscard]] LocateHit locate_walk(const Mesh& mesh, Vec3 point, std::uint32_t seed, double tol,
                                    bool require_positive_volume);

}  // namespace meshcore
