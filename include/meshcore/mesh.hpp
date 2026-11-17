#pragma once

#include "meshcore/vec3.hpp"

#include <array>
#include <cstdint>
#include <vector>

namespace meshcore {

struct Mesh {
  std::vector<Vec3> vertices;
  std::vector<std::array<std::uint32_t, 4>> tets;
  // neighbors[t][i] is the tet across the face opposite vertex i, or -1.
  std::vector<std::array<std::int32_t, 4>> neighbors;

  void build_adjacency();
  [[nodiscard]] Vec3 centroid(std::uint32_t tet) const;
};

struct Grid {
  Mesh mesh;
  int nx = 0;
  int ny = 0;
  int nz = 0;
  double h = 1.0;

  [[nodiscard]] std::uint32_t vertex_index(int i, int j, int k) const;
  [[nodiscard]] Vec3 position(int i, int j, int k) const;
  [[nodiscard]] std::uint32_t tet_index(int i, int j, int k, int local) const;
};

// Conforming 6-tet split of an nx*ny*nz cube grid. Adjacency is built.
[[nodiscard]] Grid make_cube_grid(int nx, int ny, int nz, double h = 1.0);

[[nodiscard]] Mesh make_single_tet(Vec3 a, Vec3 b, Vec3 c, Vec3 d);

// One large tet containing the sample point, and one small tet whose centroid
// is much closer to that point. Nearest-centroid search returns the small tet.
struct CentroidMiss {
  Mesh mesh;
  Vec3 point;
  std::uint32_t containing_tet = 0;
  std::uint32_t nearest_centroid_tet = 1;
};

[[nodiscard]] CentroidMiss make_centroid_miss();

[[nodiscard]] double mesh_diagonal(const Mesh& mesh);

// Barycentric weights are dimensionless, so the containment tolerance scales
// with the mesh diameter rather than with an absolute world unit.
[[nodiscard]] double absolute_tolerance(const Mesh& mesh);

}  // namespace meshcore
