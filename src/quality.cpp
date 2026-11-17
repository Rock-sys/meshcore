#include "meshcore/quality.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace meshcore {
namespace {

double face_area(Vec3 a, Vec3 b, Vec3 c) { return 0.5 * length(cross(b - a, c - a)); }

}  // namespace

TetQuality tet_quality(Vec3 v0, Vec3 v1, Vec3 v2, Vec3 v3) {
  const Vec3 e01 = v1 - v0;
  const Vec3 e02 = v2 - v0;
  const Vec3 e03 = v3 - v0;
  const Vec3 e12 = v2 - v1;
  const Vec3 e13 = v3 - v1;
  const Vec3 e23 = v3 - v2;
  const double j = scalar_triple(e01, e02, e03);

  const double len[6] = {length(e01), length(e02), length(e03), length(e12), length(e13), length(e23)};
  double min_edge = len[0];
  double max_edge = len[0];
  for (double edge : len) {
    min_edge = std::min(min_edge, edge);
    max_edge = std::max(max_edge, edge);
  }

  TetQuality q;
  q.volume = j / 6.0;
  q.edge_ratio = (min_edge == 0.0) ? std::numeric_limits<double>::infinity() : max_edge / min_edge;

  const double prod0 = len[0] * len[1] * len[2];
  const double prod1 = len[0] * len[3] * len[4];
  const double prod2 = len[1] * len[3] * len[5];
  const double prod3 = len[2] * len[4] * len[5];
  const double prod_max = std::max(std::max(prod0, prod1), std::max(prod2, prod3));
  q.scaled_jacobian = (prod_max == 0.0) ? 0.0 : std::sqrt(2.0) * j / prod_max;

  const double areas = face_area(v1, v2, v3) + face_area(v0, v2, v3) + face_area(v0, v1, v3) +
                       face_area(v0, v1, v2);
  if (j == 0.0 || areas == 0.0) {
    q.inradius = 0.0;
    q.circumradius = std::numeric_limits<double>::infinity();
    q.radius_ratio = 0.0;
    return q;
  }

  const Vec3 offset = (cross(e02, e03) * dot(e01, e01) + cross(e03, e01) * dot(e02, e02) +
                       cross(e01, e02) * dot(e03, e03)) /
                      (2.0 * j);
  q.circumradius = length(offset);
  q.inradius = 3.0 * q.volume / areas;
  q.radius_ratio = (q.circumradius == 0.0) ? 0.0 : 3.0 * q.inradius / q.circumradius;
  return q;
}

QualitySummary summarize_quality(const std::vector<Vec3>& vertices,
                                 const std::vector<std::array<std::uint32_t, 4>>& tets) {
  QualitySummary summary;
  if (tets.empty()) {
    return summary;
  }

  summary.min_edge_ratio = std::numeric_limits<double>::infinity();
  summary.min_scaled_jacobian = std::numeric_limits<double>::infinity();
  summary.min_radius_ratio = std::numeric_limits<double>::infinity();

  for (const std::array<std::uint32_t, 4>& tet : tets) {
    const TetQuality q =
        tet_quality(vertices[tet[0]], vertices[tet[1]], vertices[tet[2]], vertices[tet[3]]);
    summary.min_edge_ratio = std::min(summary.min_edge_ratio, q.edge_ratio);
    summary.max_edge_ratio = std::max(summary.max_edge_ratio, q.edge_ratio);
    summary.min_scaled_jacobian = std::min(summary.min_scaled_jacobian, q.scaled_jacobian);
    summary.min_radius_ratio = std::min(summary.min_radius_ratio, q.radius_ratio);
    if (q.volume < 0.0) {
      ++summary.inverted;
    } else if (q.volume == 0.0) {
      ++summary.degenerate;
    } else if (q.edge_ratio <= kAcceptableEdgeRatio && q.radius_ratio < kPoorRadiusRatio) {
      ++summary.slivers_missed_by_edge_ratio;
    }
  }
  return summary;
}

}  // namespace meshcore
