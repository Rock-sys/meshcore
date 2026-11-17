#pragma once

#include "meshcore/mesh.hpp"

#include <cstdint>
#include <vector>

namespace meshcore {

struct CutStats {
  std::uint64_t interior_faces = 0;
  std::uint64_t cut_faces = 0;
  std::uint32_t min_count = 0;
  std::uint32_t max_count = 0;
  double cut_ratio = 0;
  double imbalance = 0;  // max_count / (n / parts)
};

struct Partition {
  std::vector<std::uint32_t> part;
  std::uint32_t part_count = 0;
};

// part of order[i] is floor(i * parts / n). This is "cut the array into blocks".
// It follows storage order, so a renumber changes the cut.
[[nodiscard]] Partition partition_by_order(const std::vector<std::uint32_t>& order,
                                           std::uint32_t parts);

// Recursive coordinate bisection on tet centroids, splitting the longest axis
// and handing child parts in proportion to the part count. Geometry only:
// storage order does not change the cut.
[[nodiscard]] Partition partition_rcb(const Mesh& mesh, std::uint32_t parts);

// A cut face is an interior dual-graph edge whose two tets have different parts.
[[nodiscard]] CutStats measure_cut(const Mesh& mesh, const Partition& partition);

}  // namespace meshcore
