#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"
#include "../include/ndim/static_n_dim_prefix_sum.hpp"
#include <vector>
#include <stdexcept>

using namespace ndim;

TEST_CASE("StaticNDimPrefixSum 2D Operations") {
    std::vector<int> dims = {4, 4};
    StaticNDimPrefixSum<int64_t> grid(dims);

    grid.add({0, 0}, 1);
    grid.add({0, 1}, 2);
    grid.add({1, 0}, 3);
    grid.add({1, 1}, 4);

    grid.sweep();

    REQUIRE(grid.query_range({0, 0}, {0, 0}) == 1);
    REQUIRE(grid.query_range({0, 0}, {0, 1}) == 3);
    REQUIRE(grid.query_range({0, 0}, {1, 0}) == 4);
    REQUIRE(grid.query_range({0, 0}, {1, 1}) == 10);
    REQUIRE(grid.query_range({1, 1}, {1, 1}) == 4);
}

TEST_CASE("StaticNDimPrefixSum 3D Operations") {
    std::vector<int> dims = {3, 3, 3};
    StaticNDimPrefixSum<int64_t> grid(dims);

    grid.add({1, 1, 1}, 10);
    grid.add({2, 2, 2}, 20);

    grid.sweep();

    REQUIRE(grid.query_range({0, 0, 0}, {1, 1, 1}) == 10);
    REQUIRE(grid.query_range({1, 1, 1}, {2, 2, 2}) == 30);
    REQUIRE(grid.query_range({0, 0, 0}, {0, 0, 0}) == 0);
    REQUIRE(grid.query_range({2, 2, 2}, {2, 2, 2}) == 20);
}

TEST_CASE("StaticNDimPrefixSum Error Handling") {
    std::vector<int> dims = {2, 2};
    StaticNDimPrefixSum<int64_t> grid(dims);

    grid.add({0, 0}, 1);

    // Query before sweep
    REQUIRE_THROWS_AS(grid.query_range({0, 0}, {1, 1}), std::runtime_error);

    grid.sweep();

    // Add after sweep
    REQUIRE_THROWS_AS(grid.add({1, 1}, 1), std::runtime_error);
}
