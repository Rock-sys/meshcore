#pragma once

#include "meshcore/locate.hpp"

#include <cstdint>
#include <vector>

namespace meshcore {

struct Aabb {
  Vec3 lo;
  Vec3 hi;

  static Aabb from_point(Vec3 p) { return {p, p}; }
  void expand(Vec3 p);
  void expand_box(const Aabb& other);
  [[nodiscard]] bool contains(Vec3 p, double tol) const;
};

// Bounding-box tree over tetrahedra. A query walks boxes that contain the
// point, then confirms with barycentric coordinates. The mesh must outlive
// the tree and its coordinates must not change.
class AabbTree {
 public:
  explicit AabbTree(const Mesh& mesh, std::uint32_t leaf_size = 8);

  [[nodiscard]] LocateHit locate(Vec3 point, double tol, bool require_positive_volume) const;

  [[nodiscard]] std::uint32_t leaf_size() const { return leaf_size_; }

 private:
  struct Node {
    Aabb box;
    std::uint32_t left = 0;
    std::uint32_t right = 0;
    std::uint32_t begin = 0;
    std::uint32_t count = 0;
    bool leaf = false;
  };

  [[nodiscard]] std::uint32_t build(std::uint32_t begin, std::uint32_t end);
  [[nodiscard]] Aabb range_box(std::uint32_t begin, std::uint32_t end) const;
  [[nodiscard]] int split_axis(std::uint32_t begin, std::uint32_t end) const;

  const Mesh* mesh_ = nullptr;
  std::uint32_t leaf_size_ = 8;
  std::uint32_t root_ = 0;
  std::vector<std::uint32_t> order_;
  std::vector<Vec3> centroids_;
  std::vector<Node> nodes_;
};

}  // namespace meshcore
