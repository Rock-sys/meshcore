#pragma once

#include "meshcore/vec3.hpp"

#include <array>
#include <cstdint>
#include <vector>

namespace meshcore {

// Element quality for one tetrahedron.
//
// Conventions, so the numbers can be defended:
//   J = (v1-v0) · ((v2-v0) × (v3-v0)),   volume = J/6
//   edge_ratio = longest edge / shortest edge          (1 = regular)
//   scaled_jacobian = √2 * J / max_c (product of the three edge lengths at corner c)
//       1 on a regular tet, <= 0 when the tet is inverted or flat.
//       The √2 makes the regular tet score 1. Some codes omit it and top out at 1/√2.
//   radius_ratio = 3 r / R,  r = 3V / (sum of face areas),  R = circumradius
//       1 on a regular tet. Edge ratio misses slivers; radius ratio does not.
struct TetQuality {
  double volume = 0;
  double edge_ratio = 0;
  double scaled_jacobian = 0;
  double radius_ratio = 0;
  double inradius = 0;
  double circumradius = 0;
};

[[nodiscard]] TetQuality tet_quality(Vec3 v0, Vec3 v1, Vec3 v2, Vec3 v3);

// Knupp's published acceptable band for the edge-length ratio is [1, 3].
// A positive tet inside that band with a tiny radius ratio is the sliver case.
inline constexpr double kAcceptableEdgeRatio = 3.0;
inline constexpr double kPoorRadiusRatio = 0.1;

struct QualitySummary {
  double min_edge_ratio = 0;
  double max_edge_ratio = 0;
  double min_scaled_jacobian = 0;
  double min_radius_ratio = 0;
  std::uint32_t inverted = 0;
  std::uint32_t degenerate = 0;
  std::uint32_t slivers_missed_by_edge_ratio = 0;
};

[[nodiscard]] QualitySummary summarize_quality(
    const std::vector<Vec3>& vertices,
    const std::vector<std::array<std::uint32_t, 4>>& tets);

}  // namespace meshcore
