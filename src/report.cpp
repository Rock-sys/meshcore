#include "meshcore/aabb_tree.hpp"
#include "meshcore/centroid_kdtree.hpp"
#include "meshcore/locate.hpp"
#include "meshcore/mesh.hpp"
#include "meshcore/partition.hpp"
#include "meshcore/quality.hpp"
#include "meshcore/synthetic.hpp"

#include <chrono>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <vector>

namespace {

void print_quality(const char* name, const meshcore::TetQuality& q) {
  std::cout << std::left << std::setw(12) << name << "edge " << std::setw(10) << q.edge_ratio
            << " SJ " << std::setw(12) << q.scaled_jacobian << " radius " << std::setw(12)
            << q.radius_ratio << " volume " << q.volume << "\n";
}

meshcore::Vec3 interior_point(const meshcore::Mesh& mesh, std::uint32_t tet, double w0, double w1,
                             double w2, double w3) {
  const std::array<std::uint32_t, 4>& ids = mesh.tets[tet];
  return w0 * mesh.vertices[ids[0]] + w1 * mesh.vertices[ids[1]] + w2 * mesh.vertices[ids[2]] +
         w3 * mesh.vertices[ids[3]];
}

template <typename Fn>
double mean_microseconds(Fn&& fn) {
  const auto start = std::chrono::steady_clock::now();
  fn();
  const auto stop = std::chrono::steady_clock::now();
  return std::chrono::duration<double, std::micro>(stop - start).count();
}

}  // namespace

int main() {
  std::cout << std::setprecision(6);
  std::cout << "meshcore\n\n";

  meshcore::Vec3 a, b, c, d;
  meshcore::regular_tet(a, b, c, d);
  print_quality("regular", meshcore::tet_quality(a, b, c, d));
  meshcore::sliver_tet(a, b, c, d);
  print_quality("sliver", meshcore::tet_quality(a, b, c, d));
  meshcore::inverted_tet(a, b, c, d);
  print_quality("inverted", meshcore::tet_quality(a, b, c, d));
  std::cout << "sliver stays inside the edge-ratio band [1, 3] and fails radius ratio / SJ\n\n";

  const meshcore::CentroidMiss miss = meshcore::make_centroid_miss();
  const double miss_tol = meshcore::absolute_tolerance(miss.mesh);
  const meshcore::CentroidKdTree miss_tree(miss.mesh);
  const meshcore::KnnHit knn = miss_tree.knn(miss.point, 1).front();
  const auto scan = meshcore::locate_scan(miss.mesh, miss.point, miss_tol, true);
  std::cout << "centroid k=1 returns tet " << knn.tet << " (contains point: "
            << meshcore::tet_contains(miss.mesh, knn.tet, miss.point, miss_tol, true) << ")\n";
  std::cout << "scan returns tet " << *scan.tet << "\n";
  std::cout << "k=1 locate finds a tet: " << miss_tree.locate(miss.point, 1, miss_tol, true).tet.has_value()
            << "\n\n";

  constexpr int n = 16;
  const meshcore::Grid grid = meshcore::make_cube_grid(n, n, n);
  const double tol = meshcore::absolute_tolerance(grid.mesh);
  const std::uint32_t far = grid.tet_index(n - 1, n - 1, n - 1, 0);
  const meshcore::Vec3 far_point = interior_point(grid.mesh, far, 0.15, 0.2, 0.25, 0.4);
  const meshcore::LocateHit walked = meshcore::locate_walk(grid.mesh, far_point, 0, tol, true);
  const meshcore::AabbTree aabb(grid.mesh);
  const meshcore::LocateHit boxed = aabb.locate(far_point, tol, true);
  std::cout << "grid " << n << "^3  tets " << grid.mesh.tets.size() << "\n";
  std::cout << "far point walk steps " << walked.steps << ", AABB nodes " << boxed.nodes_visited
            << ", AABB barycentric tests " << boxed.barycentric_tests << "\n";

  constexpr int kSamples = 200;
  std::vector<meshcore::Vec3> points;
  points.reserve(kSamples);
  std::uint32_t rng = 7;
  for (int i = 0; i < kSamples; ++i) {
    rng = rng * 1664525u + 1013904223u;
    const std::uint32_t tet = rng % static_cast<std::uint32_t>(grid.mesh.tets.size());
    points.push_back(interior_point(grid.mesh, tet, 0.1, 0.2, 0.3, 0.4));
  }

  int mismatches = 0;
  std::uint64_t walk_steps = 0;
  std::uint64_t aabb_tests = 0;
  for (meshcore::Vec3 point : points) {
    const auto s = meshcore::locate_scan(grid.mesh, point, tol, true);
    const auto w = meshcore::locate_walk(grid.mesh, point, 0, tol, true);
    const auto x = aabb.locate(point, tol, true);
    walk_steps += w.steps;
    aabb_tests += x.barycentric_tests;
    if (!s.tet || !w.tet || !x.tet ||
        !meshcore::tet_contains(grid.mesh, *s.tet, point, tol, true) ||
        !meshcore::tet_contains(grid.mesh, *w.tet, point, tol, true) ||
        !meshcore::tet_contains(grid.mesh, *x.tet, point, tol, true)) {
      ++mismatches;
    }
  }
  std::cout << "interior samples " << kSamples << ", mismatches " << mismatches << "\n";
  std::cout << "mean walk steps " << static_cast<double>(walk_steps) / kSamples
            << ", mean AABB barycentric tests " << static_cast<double>(aabb_tests) / kSamples
            << "\n";

  for (meshcore::Vec3 point : points) {
    (void)meshcore::locate_scan(grid.mesh, point, tol, true);
    (void)meshcore::locate_walk(grid.mesh, point, 0, tol, true);
    (void)aabb.locate(point, tol, true);
  }
  const double scan_us = mean_microseconds([&] {
    for (meshcore::Vec3 point : points) {
      (void)meshcore::locate_scan(grid.mesh, point, tol, true);
    }
  });
  const double walk_us = mean_microseconds([&] {
    for (meshcore::Vec3 point : points) {
      (void)meshcore::locate_walk(grid.mesh, point, 0, tol, true);
    }
  });
  const double aabb_us = mean_microseconds([&] {
    for (meshcore::Vec3 point : points) {
      (void)aabb.locate(point, tol, true);
    }
  });
  std::cout << "mean us/query  scan " << scan_us / kSamples << "  walk " << walk_us / kSamples
            << "  AABB " << aabb_us / kSamples << "\n";

  std::cout << "AABB leaf size vs mean barycentric tests:\n";
  for (std::uint32_t leaf : {4u, 8u, 16u, 32u}) {
    const meshcore::AabbTree tree(grid.mesh, leaf);
    std::uint64_t tests = 0;
    for (meshcore::Vec3 point : points) {
      tests += tree.locate(point, tol, true).barycentric_tests;
    }
    std::cout << "  leaf " << leaf << "  tests " << static_cast<double>(tests) / kSamples << "\n";
  }
  std::cout << "\n";

  // Centroid KD-tree: index cost, query cost against a linear scan, and how
  // many candidates point location needs on this structured grid.
  const auto kd_start = std::chrono::steady_clock::now();
  const meshcore::CentroidKdTree kd(grid.mesh);
  const auto kd_stop = std::chrono::steady_clock::now();
  const double kd_build_ms = std::chrono::duration<double, std::milli>(kd_stop - kd_start).count();
  std::cout << "centroid KD-tree, " << grid.mesh.tets.size() << " tets, build " << kd_build_ms
            << " ms\n";

  for (std::uint32_t k : {1u, 8u}) {
    std::uint32_t found = 0;
    for (meshcore::Vec3 point : points) {
      if (kd.locate(point, k, tol, true).tet.has_value()) {
        ++found;
      }
    }
    std::cout << "  locate k=" << k << " found " << found << "/" << kSamples << "\n";
  }

  const double kd1_us = mean_microseconds([&] {
                          for (meshcore::Vec3 point : points) {
                            (void)kd.knn(point, 1);
                          }
                        }) / kSamples;
  const double kd8_us = mean_microseconds([&] {
                          for (meshcore::Vec3 point : points) {
                            (void)kd.knn(point, 8);
                          }
                        }) / kSamples;
  const double linear_us = mean_microseconds([&] {
                              for (meshcore::Vec3 point : points) {
                                double best = 0.0;
                                for (std::uint32_t t = 0;
                                     t < static_cast<std::uint32_t>(grid.mesh.tets.size()); ++t) {
                                  const meshcore::Vec3 delta = grid.mesh.centroid(t) - point;
                                  const double d2 = meshcore::dot(delta, delta);
                                  if (t == 0 || d2 < best) {
                                    best = d2;
                                  }
                                }
                              }
                            }) / kSamples;
  std::cout << "  mean us/query  k=1 " << kd1_us << "  k=8 " << kd8_us
            << "  single-pass linear scan " << linear_us << "  speedup " << linear_us / kd1_us
            << "x\n";

  std::cout << "  leaf size vs k=1 us/query:";
  for (std::uint32_t leaf : {4u, 8u, 16u, 32u}) {
    const meshcore::CentroidKdTree tree(grid.mesh, leaf);
    const double us = mean_microseconds([&] {
                        for (meshcore::Vec3 point : points) {
                          (void)tree.knn(point, 1);
                        }
                      }) / kSamples;
    std::cout << "  " << leaf << ":" << us;
  }
  std::cout << "\n\n";

  const meshcore::Grid part_grid = meshcore::make_cube_grid(8, 8, 8);
  const auto tn = static_cast<std::uint32_t>(part_grid.mesh.tets.size());
  std::vector<std::uint32_t> natural(tn);
  for (std::uint32_t i = 0; i < tn; ++i) {
    natural[i] = i;
  }
  std::vector<std::uint32_t> shuffled = natural;
  std::mt19937 mix(1);
  std::shuffle(shuffled.begin(), shuffled.end(), mix);
  std::cout << "dual-graph face cut, 8^3 grid, " << tn << " tets\n";
  for (std::uint32_t parts : {8u, 32u}) {
    const meshcore::CutStats natural_cut =
        meshcore::measure_cut(part_grid.mesh, meshcore::partition_by_order(natural, parts));
    const meshcore::CutStats shuffled_cut =
        meshcore::measure_cut(part_grid.mesh, meshcore::partition_by_order(shuffled, parts));
    const meshcore::CutStats rcb_cut =
        meshcore::measure_cut(part_grid.mesh, meshcore::partition_rcb(part_grid.mesh, parts));
    std::cout << "  parts " << parts << "  natural " << natural_cut.cut_ratio << "  shuffled "
              << shuffled_cut.cut_ratio << "  rcb " << rcb_cut.cut_ratio << "  rcb imbalance "
              << rcb_cut.imbalance << "\n";
  }
  return mismatches == 0 ? 0 : 1;
}
