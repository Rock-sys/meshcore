#include "meshcore/aabb_tree.hpp"

#include <algorithm>
#include <stdexcept>

namespace meshcore {

void Aabb::expand(Vec3 p) {
  lo.x = std::min(lo.x, p.x);
  lo.y = std::min(lo.y, p.y);
  lo.z = std::min(lo.z, p.z);
  hi.x = std::max(hi.x, p.x);
  hi.y = std::max(hi.y, p.y);
  hi.z = std::max(hi.z, p.z);
}

void Aabb::expand_box(const Aabb& other) {
  expand(other.lo);
  expand(other.hi);
}

bool Aabb::contains(Vec3 p, double tol) const {
  return p.x >= lo.x - tol && p.x <= hi.x + tol && p.y >= lo.y - tol && p.y <= hi.y + tol &&
         p.z >= lo.z - tol && p.z <= hi.z + tol;
}

AabbTree::AabbTree(const Mesh& mesh, std::uint32_t leaf_size) : mesh_(&mesh), leaf_size_(leaf_size) {
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
  nodes_.reserve(n * 2u + 1u);
  if (n > 0) {
    root_ = build(0, n);
  }
}

Aabb AabbTree::range_box(std::uint32_t begin, std::uint32_t end) const {
  const std::array<std::uint32_t, 4>& first = mesh_->tets[order_[begin]];
  Aabb box = Aabb::from_point(mesh_->vertices[first[0]]);
  for (std::uint32_t i = begin; i < end; ++i) {
    const std::array<std::uint32_t, 4>& tet = mesh_->tets[order_[i]];
    for (std::uint32_t id : tet) {
      box.expand(mesh_->vertices[id]);
    }
  }
  return box;
}

int AabbTree::split_axis(std::uint32_t begin, std::uint32_t end) const {
  Vec3 lo = centroids_[order_[begin]];
  Vec3 hi = lo;
  for (std::uint32_t i = begin; i < end; ++i) {
    const Vec3 c = centroids_[order_[i]];
    lo.x = std::min(lo.x, c.x);
    lo.y = std::min(lo.y, c.y);
    lo.z = std::min(lo.z, c.z);
    hi.x = std::max(hi.x, c.x);
    hi.y = std::max(hi.y, c.y);
    hi.z = std::max(hi.z, c.z);
  }
  const Vec3 span = hi - lo;
  if (span.y > span.x && span.y >= span.z) {
    return 1;
  }
  if (span.z > span.x && span.z > span.y) {
    return 2;
  }
  return 0;
}

std::uint32_t AabbTree::build(std::uint32_t begin, std::uint32_t end) {
  const std::uint32_t index = static_cast<std::uint32_t>(nodes_.size());
  nodes_.push_back(Node{});
  nodes_[index].box = range_box(begin, end);
  nodes_[index].begin = begin;
  nodes_[index].count = end - begin;

  if (end - begin <= leaf_size_) {
    nodes_[index].leaf = true;
    return index;
  }

  const int axis = split_axis(begin, end);
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
  nodes_[index].left = build(begin, mid);
  nodes_[index].right = build(mid, end);
  return index;
}

LocateHit AabbTree::locate(Vec3 point, double tol, bool require_positive_volume) const {
  LocateHit hit;
  if (nodes_.empty()) {
    return hit;
  }

  const auto visit = [&](auto&& self, std::uint32_t node) -> void {
    if (hit.tet.has_value()) {
      return;
    }
    ++hit.nodes_visited;
    const Node& current = nodes_[node];
    if (!current.box.contains(point, tol)) {
      return;
    }
    if (current.leaf) {
      for (std::uint32_t i = 0; i < current.count; ++i) {
        const std::uint32_t tet = order_[current.begin + i];
        ++hit.barycentric_tests;
        if (tet_contains(*mesh_, tet, point, tol, require_positive_volume)) {
          hit.tet = tet;
          return;
        }
      }
      return;
    }
    self(self, current.left);
    self(self, current.right);
  };
  visit(visit, root_);
  return hit;
}

}  // namespace meshcore
