#pragma once

#include "meshcore/locate.hpp"

#include <cstdint>
#include <vector>

namespace meshcore {

struct KnnHit {
  std::uint32_t tet = 0;
  double dist2 = 0;
};

// KD-tree on tet centroids. Exact for "which centroid is nearest".
// Not exact for "which tet contains this point": a large tet's centroid can
// sit far from a query that is still inside it.
class CentroidKdTree {
 public:
  explicit CentroidKdTree(const Mesh& mesh, std::uint32_t leaf_size = 8);

  [[nodiscard]] std::vector<KnnHit> knn(Vec3 point, std::uint32_t k) const;

  // Checks the k nearest centroids and returns the first tet that contains
  // the point. An empty hit means those k centroids were the wrong elements,
  // not that the point is outside the mesh.
  [[nodiscard]] LocateHit locate(Vec3 point, std::uint32_t k, double tol,
                                 bool require_positive_volume) const;

 private:
  struct Node {
    double split = 0;
    std::uint32_t left = 0;
    std::uint32_t right = 0;
    std::uint32_t begin = 0;
    std::uint32_t count = 0;
    std::uint8_t axis = 0;
    bool leaf = false;
  };

  [[nodiscard]] std::uint32_t build(std::uint32_t begin, std::uint32_t end, int depth);

  const Mesh* mesh_ = nullptr;
  std::uint32_t leaf_size_ = 8;
  std::uint32_t root_ = 0;
  std::vector<std::uint32_t> order_;
  std::vector<Vec3> centroids_;
  std::vector<Node> nodes_;
};

}  // namespace meshcore
