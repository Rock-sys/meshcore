cmake_minimum_required(VERSION 3.16)
project(meshcore LANGUAGES CXX)

# Tet quality, point location, and dual-graph RCB.
#   cmake -S . -B build
#   cmake --build build --config Release
#   ctest --test-dir build -C Release --output-on-failure
#   ./build/Release/meshcore_report   (multi-config)
#   ./build/meshcore_report           (single-config)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_library(meshcore
  src/quality.cpp
  src/mesh.cpp
  src/locate.cpp
  src/aabb_tree.cpp
  src/centroid_kdtree.cpp
  src/partition.cpp
)
target_include_directories(meshcore PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

add_executable(meshcore_tests tests/test_main.cpp)
target_link_libraries(meshcore_tests PRIVATE meshcore)

add_executable(meshcore_report src/report.cpp)
target_link_libraries(meshcore_report PRIVATE meshcore)

function(meshcore_warnings target)
  if(MSVC)
    target_compile_options(${target} PRIVATE /W4 /permissive- /utf-8)
  else()
    target_compile_options(${target} PRIVATE -Wall -Wextra -Wpedantic)
  endif()
endfunction()

meshcore_warnings(meshcore)
meshcore_warnings(meshcore_tests)
meshcore_warnings(meshcore_report)

enable_testing()
add_test(NAME meshcore_tests COMMAND $<TARGET_FILE:meshcore_tests>)
