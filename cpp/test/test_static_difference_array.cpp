#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"
#include "../include/ndim/static_n_dim_difference_array.hpp"
#include <vector>
#include <stdexcept>

using namespace ndim;

TEST_CASE("StaticNDimDifferenceArray 2D Operations") {
    std::vector<int> dims = {4, 4};
    StaticNDimDifferenceArray<int64_t> grid(dims);

    grid.add_range({0, 0}, {2, 2}, 10);
    grid.add_range({1, 1}, {3, 3}, 5);

    grid.sweep();

    REQUIRE(grid.query_range({0, 0}, {0, 0}) == 10);
    REQUIRE(grid.query_range({1, 1}, {1, 1}) == 15);
    REQUIRE(grid.query_range({3, 3}, {3, 3}) == 5);
    REQUIRE(grid.query_range({0, 0}, {3, 3}) == 10 * 9 + 5 * 9);
}

TEST_CASE("StaticNDimDifferenceArray 3D Operations") {
    std::vector<int> dims = {5, 5, 5};
    StaticNDimDifferenceArray<int64_t> grid(dims);

    grid.add_range({0, 0, 0}, {4, 4, 4}, 1);
    grid.add_range({1, 1, 1}, {3, 3, 3}, 2);

    grid.sweep();

    REQUIRE(grid.query_range({0, 0, 0}, {0, 0, 0}) == 1);
    REQUIRE(grid.query_range({2, 2, 2}, {2, 2, 2}) == 3);
    
    // (5x5x5) * 1 + (3x3x3) * 2 = 125 + 54 = 179
    REQUIRE(grid.query_range({0, 0, 0}, {4, 4, 4}) == 179);
}

TEST_CASE("StaticNDimDifferenceArray Error Handling") {
    std::vector<int> dims = {2, 2};
    StaticNDimDifferenceArray<int64_t> grid(dims);

    grid.add_range({0, 0}, {1, 1}, 1);

    // Query before sweep
    REQUIRE_THROWS_AS(grid.query_range({0, 0}, {1, 1}), std::runtime_error);

    grid.sweep();

    // Add after sweep
    REQUIRE_THROWS_AS(grid.add_range({1, 1}, {1, 1}, 1), std::runtime_error);
}
