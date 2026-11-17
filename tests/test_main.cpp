#include "meshcore/aabb_tree.hpp"
#include "meshcore/centroid_kdtree.hpp"
#include "meshcore/locate.hpp"
#include "meshcore/mesh.hpp"
#include "meshcore/partition.hpp"
#include "meshcore/quality.hpp"
#include "meshcore/synthetic.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <random>
#include <string>
#include <vector>

namespace {

int g_fails = 0;

void check(bool cond, const char* expr, const char* file, int line) {
  if (!cond) {
    std::cerr << file << ":" << line << "  " << expr << "\n";
    ++g_fails;
  }
}

void expect_close(double got, double expected, double tol, const char* name) {
  if (!(std::abs(got - expected) <= tol)) {
    std::cerr << name << " got " << got << " expected " << expected << " tol " << tol << "\n";
    ++g_fails;
  }
}

#define CHECK(expr) ::check(static_cast<bool>(expr), #expr, __FILE__, __LINE__)

meshcore::Vec3 interior_point(const meshcore::Mesh& mesh, std::uint32_t tet, double w0, double w1,
                             double w2, double w3) {
  const std::array<std::uint32_t, 4>& ids = mesh.tets[tet];
  return w0 * mesh.vertices[ids[0]] + w1 * mesh.vertices[ids[1]] + w2 * mesh.vertices[ids[2]] +
         w3 * mesh.vertices[ids[3]];
}

// Reference implementation: sort every centroid by distance.
std::vector<meshcore::KnnHit> brute_knn(const meshcore::Mesh& mesh, meshcore::Vec3 q,
                                        std::uint32_t k) {
  std::vector<meshcore::KnnHit> all;
  all.reserve(mesh.tets.size());
  for (std::uint32_t t = 0; t < static_cast<std::uint32_t>(mesh.tets.size()); ++t) {
    const meshcore::Vec3 d = mesh.centroid(t) - q;
    all.push_back(meshcore::KnnHit{t, meshcore::dot(d, d)});
  }
  std::sort(all.begin(), all.end(), [](const meshcore::KnnHit& lhs, const meshcore::KnnHit& rhs) {
    return lhs.dist2 < rhs.dist2;
  });
  if (all.size() > k) {
    all.resize(k);
  }
  return all;
}

void test_regular_metrics() {
  meshcore::Vec3 a, b, c, d;
  meshcore::regular_tet(a, b, c, d);
  const meshcore::TetQuality q = meshcore::tet_quality(a, b, c, d);
  expect_close(q.volume, 8.0 / 3.0, 1e-12, "regular volume");
  expect_close(q.edge_ratio, 1.0, 1e-12, "regular edge ratio");
  expect_close(q.scaled_jacobian, 1.0, 1e-12, "regular scaled jacobian");
  expect_close(q.radius_ratio, 1.0, 1e-12, "regular radius ratio");
  expect_close(q.circumradius, std::sqrt(3.0), 1e-12, "regular circumradius");

  const meshcore::Mesh mesh = meshcore::make_single_tet(a, b, c, d);
  const meshcore::Vec3 centroid = (a + b + c + d) * 0.25;
  const meshcore::Barycentric weights = meshcore::barycentric(a, b, c, d, centroid);
  for (double w : weights.w) {
    expect_close(w, 0.25, 1e-12, "regular centroid weight");
  }
  CHECK(meshcore::tet_contains(mesh, 0, centroid, 1e-12, true));
}

void test_sliver_hides_from_edge_ratio() {
  meshcore::Vec3 a, b, c, d;
  meshcore::sliver_tet(a, b, c, d);
  const meshcore::TetQuality q = meshcore::tet_quality(a, b, c, d);
  CHECK(q.volume > 0.0);
  CHECK(q.edge_ratio <= meshcore::kAcceptableEdgeRatio);
  CHECK(q.edge_ratio > 1.2);
  CHECK(q.scaled_jacobian > 0.0);
  CHECK(q.scaled_jacobian < 1e-3);
  CHECK(q.radius_ratio > 0.0);
  CHECK(q.radius_ratio < meshcore::kPoorRadiusRatio);

  const meshcore::Mesh mesh = meshcore::make_single_tet(a, b, c, d);
  const meshcore::QualitySummary summary = meshcore::summarize_quality(mesh.vertices, mesh.tets);
  CHECK(summary.inverted == 0);
  CHECK(summary.slivers_missed_by_edge_ratio == 1);
}

void test_edge_ratio_misses_inversion() {
  meshcore::Vec3 a, b, c, d;
  meshcore::inverted_tet(a, b, c, d);
  const meshcore::TetQuality q = meshcore::tet_quality(a, b, c, d);
  expect_close(q.edge_ratio, 1.0, 1e-12, "inverted edge ratio");
  CHECK(q.volume < 0.0);
  CHECK(q.scaled_jacobian < 0.0);
  CHECK(q.radius_ratio < 0.0);

  const meshcore::Mesh mesh = meshcore::make_single_tet(a, b, c, d);
  const meshcore::Vec3 centroid = (a + b + c + d) * 0.25;
  const double tol = meshcore::absolute_tolerance(mesh);
  CHECK(!meshcore::tet_contains(mesh, 0, centroid, tol, true));
  CHECK(meshcore::tet_contains(mesh, 0, centroid, tol, false));
  CHECK(!meshcore::locate_scan(mesh, centroid, tol, true).tet.has_value());
}

void test_grid_is_conforming() {
  const meshcore::Grid grid = meshcore::make_cube_grid(4, 4, 4);
  CHECK(grid.mesh.tets.size() == 4u * 4u * 4u * 6u);

  std::uint64_t boundary = 0;
  std::uint64_t interior = 0;
  for (std::uint32_t t = 0; t < grid.mesh.tets.size(); ++t) {
    const meshcore::TetQuality q = meshcore::tet_quality(
        grid.mesh.vertices[grid.mesh.tets[t][0]], grid.mesh.vertices[grid.mesh.tets[t][1]],
        grid.mesh.vertices[grid.mesh.tets[t][2]], grid.mesh.vertices[grid.mesh.tets[t][3]]);
    CHECK(q.volume > 0.0);
    for (std::int32_t neighbor : grid.mesh.neighbors[t]) {
      if (neighbor < 0) {
        ++boundary;
      } else if (static_cast<std::uint32_t>(neighbor) > t) {
        ++interior;
        CHECK(grid.mesh.neighbors[static_cast<std::uint32_t>(neighbor)][0] ==
                  static_cast<std::int32_t>(t) ||
              grid.mesh.neighbors[static_cast<std::uint32_t>(neighbor)][1] ==
                  static_cast<std::int32_t>(t) ||
              grid.mesh.neighbors[static_cast<std::uint32_t>(neighbor)][2] ==
                  static_cast<std::int32_t>(t) ||
              grid.mesh.neighbors[static_cast<std::uint32_t>(neighbor)][3] ==
                  static_cast<std::int32_t>(t));
      }
    }
  }
  const std::uint64_t expected_boundary = 4ull * (4ull * 4ull + 4ull * 4ull + 4ull * 4ull);
  CHECK(boundary == expected_boundary);
  CHECK(interior > 0);

  const meshcore::QualitySummary summary =
      meshcore::summarize_quality(grid.mesh.vertices, grid.mesh.tets);
  CHECK(summary.inverted == 0);
  CHECK(summary.slivers_missed_by_edge_ratio == 0);
  CHECK(summary.min_edge_ratio <= meshcore::kAcceptableEdgeRatio);
}

void test_locators_agree_inside_and_outside() {
  const meshcore::Grid grid = meshcore::make_cube_grid(6, 5, 4);
  const double tol = meshcore::absolute_tolerance(grid.mesh);
  const meshcore::AabbTree aabb(grid.mesh);
  const meshcore::CentroidKdTree kdtree(grid.mesh);

  std::uint32_t rng = 1;
  for (int sample = 0; sample < 80; ++sample) {
    rng = rng * 1664525u + 1013904223u;
    const std::uint32_t tet = rng % static_cast<std::uint32_t>(grid.mesh.tets.size());
    const meshcore::Vec3 point = interior_point(grid.mesh, tet, 0.1, 0.2, 0.3, 0.4);
    const meshcore::LocateHit scanned = meshcore::locate_scan(grid.mesh, point, tol, true);
    const meshcore::LocateHit walked = meshcore::locate_walk(grid.mesh, point, 0, tol, true);
    const meshcore::LocateHit boxed = aabb.locate(point, tol, true);
    const meshcore::LocateHit nearest = kdtree.locate(point, 32, tol, true);
    CHECK(scanned.tet.has_value());
    CHECK(walked.tet.has_value());
    CHECK(boxed.tet.has_value());
    CHECK(nearest.tet.has_value());
    CHECK(meshcore::tet_contains(grid.mesh, *scanned.tet, point, tol, true));
    CHECK(meshcore::tet_contains(grid.mesh, *walked.tet, point, tol, true));
    CHECK(meshcore::tet_contains(grid.mesh, *boxed.tet, point, tol, true));
    CHECK(*boxed.tet == tet);
    CHECK(*scanned.tet == tet);
  }

  const meshcore::Vec3 outside{-1.0, 0.5, 0.5};
  CHECK(!meshcore::locate_scan(grid.mesh, outside, tol, true).tet.has_value());
  CHECK(!meshcore::locate_walk(grid.mesh, outside, 0, tol, true).tet.has_value());
  CHECK(!aabb.locate(outside, tol, true).tet.has_value());
}

void test_centroid_knn_is_not_point_location() {
  const meshcore::CentroidMiss miss = meshcore::make_centroid_miss();
  const double tol = meshcore::absolute_tolerance(miss.mesh);
  const meshcore::CentroidKdTree kdtree(miss.mesh);
  const meshcore::AabbTree aabb(miss.mesh);

  const std::vector<meshcore::KnnHit> nearest = kdtree.knn(miss.point, 1);
  CHECK(nearest.size() == 1);
  CHECK(nearest[0].tet == miss.nearest_centroid_tet);
  CHECK(!meshcore::tet_contains(miss.mesh, nearest[0].tet, miss.point, tol, true));
  CHECK(!kdtree.locate(miss.point, 1, tol, true).tet.has_value());

  const meshcore::LocateHit recovered = kdtree.locate(miss.point, 2, tol, true);
  const meshcore::LocateHit scanned = meshcore::locate_scan(miss.mesh, miss.point, tol, true);
  const meshcore::LocateHit boxed = aabb.locate(miss.point, tol, true);
  CHECK(recovered.tet == miss.containing_tet);
  CHECK(scanned.tet == miss.containing_tet);
  CHECK(boxed.tet == miss.containing_tet);

  const meshcore::LocateHit from_container =
      meshcore::locate_walk(miss.mesh, miss.point, miss.containing_tet, tol, true);
  const meshcore::LocateHit from_distractor =
      meshcore::locate_walk(miss.mesh, miss.point, miss.nearest_centroid_tet, tol, true);
  CHECK(from_container.tet == miss.containing_tet);
  CHECK(!from_distractor.tet.has_value());
}

void test_walk_crosses_the_domain_aabb_does_not() {
  constexpr int n = 10;
  const meshcore::Grid grid = meshcore::make_cube_grid(n, n, n);
  const double tol = meshcore::absolute_tolerance(grid.mesh);
  const std::uint32_t target = grid.tet_index(n - 1, n - 1, n - 1, 0);
  const meshcore::Vec3 point = interior_point(grid.mesh, target, 0.15, 0.2, 0.25, 0.4);

  const meshcore::LocateHit walked = meshcore::locate_walk(grid.mesh, point, 0, tol, true);
  const meshcore::AabbTree aabb(grid.mesh);
  const meshcore::LocateHit boxed = aabb.locate(point, tol, true);
  CHECK(walked.tet == target);
  CHECK(boxed.tet == target);
  CHECK(walked.barycentric_tests > boxed.barycentric_tests);
  CHECK(walked.steps >= static_cast<std::uint32_t>(n));
  std::cout << "  far-point walk_steps=" << walked.steps
            << " aabb_nodes=" << boxed.nodes_visited
            << " aabb_tests=" << boxed.barycentric_tests << "\n";
}

void test_rcb_ignores_storage_order() {
  const meshcore::Grid grid = meshcore::make_cube_grid(8, 8, 8);
  const auto n = static_cast<std::uint32_t>(grid.mesh.tets.size());
  std::vector<std::uint32_t> natural(n);
  for (std::uint32_t i = 0; i < n; ++i) {
    natural[i] = i;
  }
  std::vector<std::uint32_t> shuffled = natural;
  std::mt19937 rng(1);
  std::shuffle(shuffled.begin(), shuffled.end(), rng);

  for (std::uint32_t parts : {8u, 32u}) {
    const meshcore::Partition linear = meshcore::partition_by_order(natural, parts);
    const meshcore::Partition mixed = meshcore::partition_by_order(shuffled, parts);
    const meshcore::Partition rcb = meshcore::partition_rcb(grid.mesh, parts);
    const meshcore::Partition rcb_again = meshcore::partition_rcb(grid.mesh, parts);
    CHECK(rcb.part == rcb_again.part);

    const meshcore::CutStats natural_cut = meshcore::measure_cut(grid.mesh, linear);
    const meshcore::CutStats shuffled_cut = meshcore::measure_cut(grid.mesh, mixed);
    const meshcore::CutStats rcb_cut = meshcore::measure_cut(grid.mesh, rcb);
    CHECK(rcb_cut.max_count - rcb_cut.min_count <= 1);
    CHECK(rcb_cut.min_count >= 1);
    CHECK(shuffled_cut.cut_ratio > 0.7);
    CHECK(rcb_cut.cut_ratio < (parts == 8 ? 0.10 : 0.20));
    CHECK(shuffled_cut.cut_faces > rcb_cut.cut_faces);
    CHECK(shuffled_cut.cut_ratio > natural_cut.cut_ratio);
    std::cout << "  parts=" << parts << " natural_cut=" << natural_cut.cut_ratio
              << " shuffled_cut=" << shuffled_cut.cut_ratio << " rcb_cut=" << rcb_cut.cut_ratio
              << " rcb_imbalance=" << rcb_cut.imbalance << "\n";
  }
}

// Pruning is the only part of a kd-tree that can silently return a mostly
// right answer: a bad bound drops the true nearest centroid and no other test
// in the suite would notice.
void test_knn_matches_brute_force() {
  const meshcore::Grid grid = meshcore::make_cube_grid(10, 10, 10);
  const meshcore::CentroidKdTree kdtree(grid.mesh);
  const auto n = static_cast<std::uint32_t>(grid.mesh.tets.size());

  std::uint32_t rng = 2024;
  for (int sample = 0; sample < 400; ++sample) {
    rng = rng * 1664525u + 1013904223u;
    const std::uint32_t tet = rng % n;
    // Weights must sum to one; otherwise the point is not inside the tet.
    const double w0 = 1.0 + static_cast<double>((rng >> 3) % 97);
    const double w1 = 1.0 + static_cast<double>((rng >> 9) % 97);
    const double w2 = 1.0 + static_cast<double>((rng >> 15) % 97);
    const double w3 = 1.0 + static_cast<double>((rng >> 21) % 97);
    const double sum = w0 + w1 + w2 + w3;
    const meshcore::Vec3 point =
        interior_point(grid.mesh, tet, w0 / sum, w1 / sum, w2 / sum, w3 / sum);

    for (std::uint32_t k : {1u, 8u, 32u}) {
      const std::vector<meshcore::KnnHit> got = kdtree.knn(point, k);
      const std::vector<meshcore::KnnHit> want = brute_knn(grid.mesh, point, k);
      CHECK(got.size() == want.size());
      std::vector<double> got_dist;
      std::vector<double> want_dist;
      for (const meshcore::KnnHit& hit : got) {
        got_dist.push_back(hit.dist2);
      }
      for (const meshcore::KnnHit& hit : want) {
        want_dist.push_back(hit.dist2);
      }
      std::sort(got_dist.begin(), got_dist.end());
      std::sort(want_dist.begin(), want_dist.end());
      for (std::size_t i = 0; i < std::min(got_dist.size(), want_dist.size()); ++i) {
        expect_close(got_dist[i], want_dist[i], 1e-9, "knn neighbour distance");
      }
      for (std::size_t i = 0; i < got.size(); ++i) {
        for (std::size_t j = i + 1; j < got.size(); ++j) {
          CHECK(got[i].tet != got[j].tet);
        }
      }
    }
  }
}

void test_knn_edge_cases() {
  const meshcore::Mesh empty;
  const meshcore::CentroidKdTree empty_tree(empty);
  CHECK(empty_tree.knn({0.0, 0.0, 0.0}, 1).empty());
  CHECK(!empty_tree.locate({0.0, 0.0, 0.0}, 1, 1e-12, true).tet.has_value());

  meshcore::Vec3 a, b, c, d;
  meshcore::regular_tet(a, b, c, d);
  const meshcore::Mesh one = meshcore::make_single_tet(a, b, c, d);
  const meshcore::CentroidKdTree one_tree(one);
  // Asking for more neighbours than exist must clamp rather than throw.
  CHECK(one_tree.knn(meshcore::Vec3{0.0, 0.0, 0.0}, 5).size() == 1);
}

void run(const char* name, void (*fn)()) {
  const int before = g_fails;
  std::cout << name << "\n";
  fn();
  if (g_fails == before) {
    std::cout << "  PASS\n";
  } else {
    std::cout << "  FAIL\n";
  }
}

}  // namespace

int main() {
  run("regular tet scores 1 / 1 / 1", test_regular_metrics);
  run("edge ratio accepts a sliver", test_sliver_hides_from_edge_ratio);
  run("edge ratio accepts an inversion", test_edge_ratio_misses_inversion);
  run("cube grid is conforming and positive", test_grid_is_conforming);
  run("scan, walk, and AABB agree", test_locators_agree_inside_and_outside);
  run("nearest centroid is not the containing tet", test_centroid_knn_is_not_point_location);
  run("walk crosses the mesh, AABB does not", test_walk_crosses_the_domain_aabb_does_not);
  run("RCB cut survives a renumber", test_rcb_ignores_storage_order);
  run("knn matches brute force on random queries", test_knn_matches_brute_force);
  run("knn clamps k and survives an empty mesh", test_knn_edge_cases);
  if (g_fails != 0) {
    std::cerr << g_fails << " check(s) failed\n";
    return 1;
  }
  std::cout << "all checks passed\n";
  return 0;
}
