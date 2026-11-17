#include "meshcore/mesh.hpp"

#include <algorithm>
#include <array>
#include <cstdint>
#include <stdexcept>
#include <unordered_map>
#include <utility>

namespace meshcore {
namespace {

struct FaceKey {
  std::uint32_t a = 0;
  std::uint32_t b = 0;
  std::uint32_t c = 0;

  friend bool operator==(FaceKey lhs, FaceKey rhs) {
    return lhs.a == rhs.a && lhs.b == rhs.b && lhs.c == rhs.c;
  }
};

struct FaceHash {
  std::size_t operator()(FaceKey key) const {
    std::size_t h = key.a;
    h ^= key.b + 0x9e3779b9u + (h << 6) + (h >> 2);
    h ^= key.c + 0x9e3779b9u + (h << 6) + (h >> 2);
    return h;
  }
};

FaceKey make_key(std::uint32_t x, std::uint32_t y, std::uint32_t z) {
  if (x > y) {
    std::swap(x, y);
  }
  if (y > z) {
    std::swap(y, z);
  }
  if (x > y) {
    std::swap(x, y);
  }
  return {x, y, z};
}

constexpr int kOppositeFace[4][3] = {{1, 2, 3}, {0, 2, 3}, {0, 1, 3}, {0, 1, 2}};

}  // namespace

void Mesh::build_adjacency() {
  neighbors.assign(tets.size(), std::array<std::int32_t, 4>{-1, -1, -1, -1});

  struct Hit {
    std::uint32_t tet = 0;
    std::uint8_t face = 0;
  };
  std::unordered_map<FaceKey, std::vector<Hit>, FaceHash> faces;
  faces.reserve(tets.size() * 2);

  for (std::uint32_t t = 0; t < static_cast<std::uint32_t>(tets.size()); ++t) {
    for (std::uint8_t slot = 0; slot < 4; ++slot) {
      const FaceKey key = make_key(tets[t][kOppositeFace[slot][0]], tets[t][kOppositeFace[slot][1]],
                                   tets[t][kOppositeFace[slot][2]]);
      std::vector<Hit>& hits = faces[key];
      if (hits.size() == 2) {
        throw std::runtime_error("non-manifold face: more than two tets share a triangle");
      }
      hits.push_back(Hit{t, slot});
    }
  }

  for (const auto& entry : faces) {
    const std::vector<Hit>& hits = entry.second;
    if (hits.size() == 2) {
      neighbors[hits[0].tet][hits[0].face] = static_cast<std::int32_t>(hits[1].tet);
      neighbors[hits[1].tet][hits[1].face] = static_cast<std::int32_t>(hits[0].tet);
    }
  }
}

Vec3 Mesh::centroid(std::uint32_t tet) const {
  const std::array<std::uint32_t, 4>& ids = tets[tet];
  return (vertices[ids[0]] + vertices[ids[1]] + vertices[ids[2]] + vertices[ids[3]]) * 0.25;
}

std::uint32_t Grid::vertex_index(int i, int j, int k) const {
  return (static_cast<std::uint32_t>(k) * static_cast<std::uint32_t>(ny + 1) +
          static_cast<std::uint32_t>(j)) *
             static_cast<std::uint32_t>(nx + 1) +
         static_cast<std::uint32_t>(i);
}

Vec3 Grid::position(int i, int j, int k) const {
  return {i * h, j * h, k * h};
}

std::uint32_t Grid::tet_index(int i, int j, int k, int local) const {
  const auto cell = (static_cast<std::uint32_t>(k) * static_cast<std::uint32_t>(ny) +
                     static_cast<std::uint32_t>(j)) *
                        static_cast<std::uint32_t>(nx) +
                    static_cast<std::uint32_t>(i);
  return cell * 6u + static_cast<std::uint32_t>(local);
}

Grid make_cube_grid(int nx, int ny, int nz, double h) {
  if (nx <= 0 || ny <= 0 || nz <= 0 || h <= 0.0) {
    throw std::invalid_argument("grid dimensions must be positive");
  }

  Grid grid;
  grid.nx = nx;
  grid.ny = ny;
  grid.nz = nz;
  grid.h = h;

  const int vx = nx + 1;
  const int vy = ny + 1;
  const int vz = nz + 1;
  grid.mesh.vertices.reserve(static_cast<std::size_t>(vx) * static_cast<std::size_t>(vy) *
                             static_cast<std::size_t>(vz));
  for (int k = 0; k < vz; ++k) {
    for (int j = 0; j < vy; ++j) {
      for (int i = 0; i < vx; ++i) {
        grid.mesh.vertices.push_back(grid.position(i, j, k));
      }
    }
  }

  // Freudenthal split. The shared diagonal on every face of a cube is the one
  // that joins the lexicographically low corner to the high corner, so adjacent
  // cubes agree on the diagonal and the mesh is conforming.
  constexpr int kLocal[6][4] = {
      {0, 1, 3, 7}, {0, 1, 7, 5}, {0, 3, 2, 7}, {0, 2, 6, 7}, {0, 6, 4, 7}, {0, 4, 5, 7},
  };
  constexpr int kCorner[8][3] = {
      {0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {1, 1, 0}, {0, 0, 1}, {1, 0, 1}, {0, 1, 1}, {1, 1, 1},
  };

  grid.mesh.tets.reserve(static_cast<std::size_t>(nx) * static_cast<std::size_t>(ny) *
                         static_cast<std::size_t>(nz) * 6u);
  for (int k = 0; k < nz; ++k) {
    for (int j = 0; j < ny; ++j) {
      for (int i = 0; i < nx; ++i) {
        std::uint32_t corner[8];
        for (int c = 0; c < 8; ++c) {
          corner[c] = grid.vertex_index(i + kCorner[c][0], j + kCorner[c][1], k + kCorner[c][2]);
        }
        for (const auto& local : kLocal) {
          const std::array<std::uint32_t, 4> tet = {corner[local[0]], corner[local[1]],
                                                    corner[local[2]], corner[local[3]]};
          const double jvol = scalar_triple(grid.mesh.vertices[tet[1]] - grid.mesh.vertices[tet[0]],
                                            grid.mesh.vertices[tet[2]] - grid.mesh.vertices[tet[0]],
                                            grid.mesh.vertices[tet[3]] - grid.mesh.vertices[tet[0]]);
          if (jvol <= 0.0) {
            throw std::runtime_error("cube split produced a non-positive tet");
          }
          grid.mesh.tets.push_back(tet);
        }
      }
    }
  }

  grid.mesh.build_adjacency();
  return grid;
}

Mesh make_single_tet(Vec3 a, Vec3 b, Vec3 c, Vec3 d) {
  Mesh mesh;
  mesh.vertices = {a, b, c, d};
  mesh.tets.push_back({0, 1, 2, 3});
  mesh.build_adjacency();
  return mesh;
}

CentroidMiss make_centroid_miss() {
  CentroidMiss miss;
  miss.point = {1.0, 1.0, 1.0};
  miss.containing_tet = 0;
  miss.nearest_centroid_tet = 1;

  miss.mesh.vertices = {
      {0.0, 0.0, 0.0},
      {20.0, 0.0, 0.0},
      {0.0, 20.0, 0.0},
      {0.0, 0.0, 20.0},
      {1.2, 1.2, 1.2},
      {1.5, 1.2, 1.2},
      {1.2, 1.5, 1.2},
      {1.2, 1.2, 1.5},
  };
  miss.mesh.tets.push_back({0, 1, 2, 3});
  miss.mesh.tets.push_back({4, 5, 6, 7});
  miss.mesh.build_adjacency();
  return miss;
}

double mesh_diagonal(const Mesh& mesh) {
  if (mesh.vertices.empty()) {
    return 0.0;
  }
  Vec3 lo = mesh.vertices.front();
  Vec3 hi = lo;
  for (Vec3 v : mesh.vertices) {
    lo.x = std::min(lo.x, v.x);
    lo.y = std::min(lo.y, v.y);
    lo.z = std::min(lo.z, v.z);
    hi.x = std::max(hi.x, v.x);
    hi.y = std::max(hi.y, v.y);
    hi.z = std::max(hi.z, v.z);
  }
  return length(hi - lo);
}

double absolute_tolerance(const Mesh& mesh) {
  return std::max(1.0e-14, mesh_diagonal(mesh) * 1.0e-12);
}

}  // namespace meshcore
