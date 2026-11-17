#include "meshcore/centroid_kdtree.hpp"

#include <algorithm>
#include <queue>
#include <stdexcept>

namespace meshcore {

CentroidKdTree::CentroidKdTree(const Mesh& mesh, std::uint32_t leaf_size)
    : mesh_(&mesh), leaf_size_(leaf_size) {
  if (leaf_size_ == 0) {
    throw std::invalid_argument("leaf size must be positive");
  }
  const auto n = static_cast<std::uint32_t>(mesh.tets.size());
  order_.resize(n);
  centroids_.resize(n);
  for (std::uint32_t t = 0; t < n; ++t) {
    order_[t] = t;
    centroids_[t] = mesh.centroid(t);
  }
  // A balanced tree over n items with at most leaf_size per leaf has about
  // 2*ceil(n/leaf)-1 nodes. Reserving this instead of n*2+1 keeps the index
  // proportional to the input on large meshes (~197 KB vs ~1.5 MB at 24k tets).
  const std::uint32_t leaves = (n + leaf_size_ - 1u) / leaf_size_;
  nodes_.reserve(2u * leaves + 2u);
  if (n > 0) {
    root_ = build(0, n, 0);
  }
}

std::uint32_t CentroidKdTree::build(std::uint32_t begin, std::uint32_t end, int depth) {
  const std::uint32_t index = static_cast<std::uint32_t>(nodes_.size());
  nodes_.push_back(Node{});
  nodes_[index].begin = begin;
  nodes_[index].count = end - begin;
  nodes_[index].axis = static_cast<std::uint8_t>(depth % 3);

  if (end - begin <= leaf_size_) {
    nodes_[index].leaf = true;
    return index;
  }

  const int axis = depth % 3;
  const std::uint32_t mid = begin + (end - begin) / 2u;
  std::nth_element(order_.begin() + begin, order_.begin() + mid, order_.begin() + end,
                   [&](std::uint32_t lhs, std::uint32_t rhs) {
                     const double lc = coord(centroids_[lhs], axis);
                     const double rc = coord(centroids_[rhs], axis);
                     if (lc < rc) {
                       return true;
                     }
                     if (rc < lc) {
                       return false;
                     }
                     return lhs < rhs;
                   });
  nodes_[index].split = coord(centroids_[order_[mid]], axis);
  // Both recursive calls may grow nodes_ and invalidate any reference taken
  // before them; the assigns therefore come after, re-indexing from scratch.
  const std::uint32_t left = build(begin, mid, depth + 1);
  const std::uint32_t right = build(mid, end, depth + 1);
  nodes_[index].left = left;
  nodes_[index].right = right;
  return index;
}

std::vector<KnnHit> CentroidKdTree::knn(Vec3 point, std::uint32_t k) const {
  if (k == 0) {
    throw std::invalid_argument("k must be positive");
  }
  std::vector<KnnHit> result;
  if (nodes_.empty()) {
    return result;
  }
  k = std::min(k, static_cast<std::uint32_t>(order_.size()));

  struct Item {
    double dist2 = 0;
    std::uint32_t tet = 0;
  };
  struct FarthestFirst {
    bool operator()(const Item& lhs, const Item& rhs) const { return lhs.dist2 < rhs.dist2; }
  };
  std::priority_queue<Item, std::vector<Item>, FarthestFirst> worst;

  const auto consider = [&](std::uint32_t tet) {
    const double dist2 = dot(centroids_[tet] - point, centroids_[tet] - point);
    if (worst.size() == k && dist2 >= worst.top().dist2) {
      return;
    }
    if (worst.size() == k) {
      worst.pop();
    }
    worst.push(Item{dist2, tet});
  };

  const auto visit = [&](auto&& self, std::uint32_t node) -> void {
    const Node& current = nodes_[node];
    if (current.leaf) {
      for (std::uint32_t i = 0; i < current.count; ++i) {
        consider(order_[current.begin + i]);
      }
      return;
    }
    const double q = coord(point, current.axis);
    const bool go_left = q < current.split;
    const std::uint32_t near = go_left ? current.left : current.right;
    const std::uint32_t far = go_left ? current.right : current.left;
    self(self, near);
    const double plane = q - current.split;
    const bool want_far = worst.size() < k || plane * plane < worst.top().dist2;
    if (want_far) {
      self(self, far);
    }
  };
  visit(visit, root_);

  result.resize(worst.size());
  for (std::size_t i = result.size(); i-- > 0;) {
    result[i] = KnnHit{worst.top().tet, worst.top().dist2};
    worst.pop();
  }
  return result;
}

LocateHit CentroidKdTree::locate(Vec3 point, std::uint32_t k, double tol,
                                 bool require_positive_volume) const {
  LocateHit hit;
  const std::vector<KnnHit> nearest = knn(point, k);
  for (const KnnHit& candidate : nearest) {
    ++hit.barycentric_tests;
    if (tet_contains(*mesh_, candidate.tet, point, tol, require_positive_volume)) {
      hit.tet = candidate.tet;
      return hit;
    }
  }
  return hit;
}

}  // namespace meshcore
