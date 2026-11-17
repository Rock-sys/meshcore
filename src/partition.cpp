#include "meshcore/partition.hpp"

#include <algorithm>
#include <stdexcept>

namespace meshcore {
namespace {

int longest_axis(const std::vector<Vec3>& centroids, const std::vector<std::uint32_t>& order,
                 std::uint32_t begin, std::uint32_t end) {
  Vec3 lo = centroids[order[begin]];
  Vec3 hi = lo;
  for (std::uint32_t i = begin; i < end; ++i) {
    const Vec3 c = centroids[order[i]];
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

}  // namespace

Partition partition_by_order(const std::vector<std::uint32_t>& order, std::uint32_t parts) {
  const auto n = static_cast<std::uint32_t>(order.size());
  if (parts == 0 || n == 0 || parts > n) {
    throw std::invalid_argument("parts must be in [1, n]");
  }
  Partition out;
  out.part.assign(n, 0);
  out.part_count = parts;
  for (std::uint32_t i = 0; i < n; ++i) {
    const std::uint32_t tet = order[i];
    if (tet >= n) {
      throw std::invalid_argument("order entry is out of range");
    }
    out.part[tet] = static_cast<std::uint32_t>((static_cast<std::uint64_t>(i) * parts) / n);
  }
  return out;
}

Partition partition_rcb(const Mesh& mesh, std::uint32_t parts) {
  const auto n = static_cast<std::uint32_t>(mesh.tets.size());
  if (parts == 0 || n == 0 || parts > n) {
    throw std::invalid_argument("parts must be in [1, n]");
  }

  std::vector<Vec3> centroids(n);
  std::vector<std::uint32_t> order(n);
  for (std::uint32_t t = 0; t < n; ++t) {
    centroids[t] = mesh.centroid(t);
    order[t] = t;
  }

  Partition out;
  out.part.assign(n, 0);
  out.part_count = parts;

  struct Builder {
    const std::vector<Vec3>& centroids;
    std::vector<std::uint32_t>& order;
    std::vector<std::uint32_t>& part;
    std::uint32_t next = 0;

    void run(std::uint32_t begin, std::uint32_t end, std::uint32_t k) {
      if (k == 1 || end - begin == 1) {
        const std::uint32_t id = next++;
        for (std::uint32_t i = begin; i < end; ++i) {
          part[order[i]] = id;
        }
        return;
      }
      const std::uint32_t left_k = k / 2u;
      const std::uint32_t right_k = k - left_k;
      std::uint32_t left_n =
          static_cast<std::uint32_t>((static_cast<std::uint64_t>(end - begin) * left_k) / k);
      if (left_n == 0) {
        left_n = 1;
      }
      if (left_n >= end - begin) {
        left_n = end - begin - 1;
      }
      const int axis = longest_axis(centroids, order, begin, end);
      const std::uint32_t mid = begin + left_n;
      std::nth_element(order.begin() + begin, order.begin() + mid, order.begin() + end,
                       [&](std::uint32_t lhs, std::uint32_t rhs) {
                         const double lc = coord(centroids[lhs], axis);
                         const double rc = coord(centroids[rhs], axis);
                         if (lc < rc) {
                           return true;
                         }
                         if (rc < lc) {
                           return false;
                         }
                         return lhs < rhs;
                       });
      run(begin, mid, left_k);
      run(mid, end, right_k);
    }
  };

  Builder builder{centroids, order, out.part, 0};
  builder.run(0, n, parts);
  return out;
}

CutStats measure_cut(const Mesh& mesh, const Partition& partition) {
  if (mesh.neighbors.size() != mesh.tets.size()) {
    throw std::logic_error("measure_cut requires adjacency");
  }
  if (partition.part.size() != mesh.tets.size() || partition.part_count == 0) {
    throw std::invalid_argument("partition does not match the mesh");
  }

  CutStats stats;
  std::vector<std::uint32_t> counts(partition.part_count, 0);
  for (std::uint32_t t = 0; t < partition.part.size(); ++t) {
    const std::uint32_t part = partition.part[t];
    if (part >= partition.part_count) {
      throw std::invalid_argument("part id is out of range");
    }
    ++counts[part];
    for (std::int32_t neighbor : mesh.neighbors[t]) {
      if (neighbor > static_cast<std::int32_t>(t)) {
        ++stats.interior_faces;
        if (partition.part[static_cast<std::uint32_t>(neighbor)] != part) {
          ++stats.cut_faces;
        }
      }
    }
  }

  stats.min_count = *std::min_element(counts.begin(), counts.end());
  stats.max_count = *std::max_element(counts.begin(), counts.end());
  const double average =
      static_cast<double>(mesh.tets.size()) / static_cast<double>(partition.part_count);
  stats.imbalance = stats.max_count / average;
  stats.cut_ratio =
      stats.interior_faces == 0
          ? 0.0
          : static_cast<double>(stats.cut_faces) / static_cast<double>(stats.interior_faces);
  return stats;
}

}  // namespace meshcore
