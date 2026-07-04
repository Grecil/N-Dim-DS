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

    REQUIRE_THROWS_AS(grid.query_range({0, 0}, {1, 1}), std::runtime_error);

    grid.sweep();

    REQUIRE_THROWS_AS(grid.add({1, 1}, 1), std::runtime_error);
}

#include <random>
#include <map>

TEST_CASE("StaticNDimPrefixSum Stress Test") {
    vector<int> dims = {4, 4, 4, 4, 4};
    StaticNDimPrefixSum<int64_t> grid(dims);
    map<vector<int>, int64_t> grid_bf;

    mt19937 rng(42);
    uniform_int_distribution<int> dist_c(0, 3);
    uniform_int_distribution<int64_t> dist_v(1, 100);

    for (int i = 0; i < 1000; ++i) {
        vector<int> coords = {dist_c(rng), dist_c(rng), dist_c(rng), dist_c(rng), dist_c(rng)};
        int64_t val = dist_v(rng);
        grid.add(coords, val);
        grid_bf[coords] += val;
    }

    grid.sweep();

    for (int i = 0; i < 1000; ++i) {
        vector<int> x_coords(5), y_coords(5);
        for (int d = 0; d < 5; ++d) {
            int a = dist_c(rng);
            int b = dist_c(rng);
            x_coords[d] = min(a, b);
            y_coords[d] = max(a, b);
        }

        int64_t expected = 0;
        for (int c0 = x_coords[0]; c0 <= y_coords[0]; ++c0) {
            for (int c1 = x_coords[1]; c1 <= y_coords[1]; ++c1) {
                for (int c2 = x_coords[2]; c2 <= y_coords[2]; ++c2) {
                    for (int c3 = x_coords[3]; c3 <= y_coords[3]; ++c3) {
                        for (int c4 = x_coords[4]; c4 <= y_coords[4]; ++c4) {
                            expected += grid_bf[{c0, c1, c2, c3, c4}];
                        }
                    }
                }
            }
        }
        REQUIRE(grid.query_range(x_coords, y_coords) == expected);
    }
}
