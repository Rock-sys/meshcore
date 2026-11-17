# meshcore

Small, self-contained C++17 toolkit for the three things every tetrahedral
simulation code needs before it can trust a mesh: **element quality**, **point
location**, and **partitioning**. Written as groundwork for 3D solver work
(finite-element / finite-volume, mesh generators) rather than as a general
meshing library.

Everything here is deterministic, header-plus-source, and has no dependencies
beyond the standard library.

## Modules

| Module | What it does |
| --- | --- |
| `quality` | Volume, edge ratio, scaled Jacobian, radius ratio, inradius / circumradius for one tet, plus a mesh-level summary that counts inverted and degenerate elements. |
| `mesh` | Conforming six-tet split of a structured cube grid (the Freudenthal / Kuhn subdivision), plus a face-adjacency map over shared faces. |
| `locate` | Barycentric coordinates, a straight walk across the face of the most negative weight, and a brute-force scan as ground truth. |
| `aabb_tree` | Bounding-box tree over tet boxes; walks boxes that contain the query, then confirms with a barycentric test. |
| `centroid_kdtree` | KD-tree over tet centroids: exact k-nearest-element queries, and a cheap candidate-then-confirm locator. |
| `partition` | Recursive coordinate bisection over tet centroids, dual-graph cut measurement, and a storage-order baseline to compare against. |

## Build

Windows (MSVC), from the project root:

```bat
cmake -S . -B build
cmake --build build --config Release
ctest --test-dir build -C Release --output-on-failure
build\Release\meshcore_report.exe
```

Linux / macOS (GCC or Clang):

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build --output-on-failure
./build/meshcore_report
```

C++17, no third-party dependencies. Warnings are on at `/W4 /permissive-` for
MSVC and `-Wall -Wextra -Wpedantic` elsewhere.

## Tests

`ctest` runs `meshcore_tests`, ten cases that check numbers rather than
behaviour descriptions:

```
regular tet scores 1 / 1 / 1
edge ratio accepts a sliver
edge ratio accepts an inversion
cube grid is conforming and positive
scan, walk, and AABB agree
nearest centroid is not the containing tet
walk crosses the mesh, AABB does not
RCB cut survives a renumber
knn matches brute force on random queries
knn clamps k and survives an empty mesh
```

The last two are property tests: 400 random interior query points at k = 1, 8,
32 are compared against a brute-force reference, and the returned elements must
be distinct.

## Reference numbers

`meshcore_report` on a 16 x 16 x 16 cube grid (24,576 tets), Release build,
single thread, mean over 200 interior query points:

| Locator | Build | us / query | Mean work |
| --- | --- | --- | --- |
| Brute-force scan | none | 132.5 | 24,576 barycentric tests |
| Walk | none | 1.70 | 39.4 steps |
| AABB tree | ~0.4 ms | 0.40 | 3.57 barycentric tests |
| Centroid KD-tree (k = 1) | 5.3 ms | 0.28 | 1 test |

All four agree on the 200 interior samples: zero mismatches.

Nearest-centroid queries, k = 1, against a single-pass linear scan on the same
queries:

```
centroid KD-tree, 24576 tets, build 5.2666 ms
  locate k=1 found 200/200
  locate k=8 found 200/200
  mean us/query  k=1 0.2775  k=8 1.0445  single-pass linear scan 64.801  speedup 233.517x
  leaf size vs k=1 us/query:  4:0.376  8:0.284  16:0.4155  32:1.045
```

Leaf size 8 is the default because it wins on both trees: 0.28 us per k = 1
query, and the fewest barycentric tests per AABB locate.

Absolute numbers are machine-dependent; the ratios are the point, and they are
reproduced by `meshcore_report`, which needs no external data.

## Conventions worth defending

- `scaled_jacobian = sqrt(2) * J / max_c(product of the three edge lengths at
  corner c)`. The `sqrt(2)` puts a regular tet at exactly 1.0. Codes that omit
  it top out at 1/sqrt(2), so their "0.7 is fine" band means something else.
- Edge ratio cannot see slivers: the equilateral-base, near-flat-apex tet in
  `synthetic.hpp` has edge ratio below 2 and a radius ratio of 2.1e-4.
  `summarize_quality` reports both `slivers_missed_by_edge_ratio` and the
  inverted count, because those are the two failures a solver actually pays for.
- Containment tolerance scales with the mesh diagonal, not with world units:
  crossing a 10 m mesh at `1e-12` absolute is not the same test as crossing a
  10 um mesh at `1e-12`.
- Location and orientation are separate predicates. An inverted tet still has a
  geometric interior, so `require_positive_volume` is opt-in.
- The RCB cut is a function of geometry only. The test renumbers tets and
  requires the cut to be identical; the storage-order baseline is included to
  show how badly it fails under the same renumber (0.156 -> 0.874 cut ratio).

## Known limitations

These are real and are why the numbers above carry a sample size.

- `CentroidKdTree::locate` is **not** an exact point locator. A large element's
  centroid can sit far from a point that is still inside it, so a k = 1 search
  can miss. On 2,000 random queries over the same grid it found the containing
  tet 93.5% of the time for k = 1 and 100% for k = 8. That case is preserved as
  the test `nearest centroid is not the containing tet`.
- The walk needs a seed inside the same connected component and a hop limit;
  a sliver makes the orientation test unreliable and the walk can cycle.
- Input is generated grids only. There is no `.msh` / `.vtk` reader, no
  parallel build, and no memory budget tuning beyond reserving node capacity
  from the expected leaf count.
- Quality metrics are correct but not vectorised; the KD-tree is the only
  structure with a measured query budget.

## Layout

```
include/meshcore/   headers (one per module)
src/                implementations, plus report.cpp with the benchmarks
tests/test_main.cpp ten cases, no test framework dependency
CMakeLists.txt      C++17, warnings on, ctest wired up
```

## License

MIT. See `LICENSE`.