#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"
#include "../include/ndim/dynamic_n_dim_seg_tree.hpp"
#include <vector>
#include <algorithm>
#include <limits>
#include <numeric>

using namespace ndim;

TEST_CASE("DynamicNDimSegTree 2D Operations Min") {
    auto min_func = [](int64_t a, int64_t b) { return std::min(a, b); };
    int64_t def = std::numeric_limits<int64_t>::max();
    DynamicNDimSegTree<int64_t, decltype(min_func)> tree({5, 5}, min_func, def);
    
    tree.update({1, 1}, 10);
    tree.update({2, 2}, 5);
    tree.update({3, 3}, 20);

    REQUIRE(tree.query_range({0, 0}, {4, 4}) == 5);
    REQUIRE(tree.query_range({0, 0}, {1, 1}) == 10);
    REQUIRE(tree.query_range({3, 3}, {4, 4}) == 20);
}

TEST_CASE("DynamicNDimSegTree 3D Operations GCD") {
    auto gcd_func = [](int64_t a, int64_t b) { return std::gcd(a, b); };
    int64_t def = 0;
    DynamicNDimSegTree<int64_t, decltype(gcd_func)> tree({4, 4, 4}, gcd_func, def);
    
    tree.update({0, 0, 0}, 12);
    tree.update({1, 1, 1}, 18);
    tree.update({2, 2, 2}, 24);

    REQUIRE(tree.query_range({0, 0, 0}, {1, 1, 1}) == 6);
    REQUIRE(tree.query_range({1, 1, 1}, {2, 2, 2}) == 6);
    REQUIRE(tree.query_range({0, 0, 0}, {2, 2, 2}) == 6);
    REQUIRE(tree.query_range({2, 2, 2}, {3, 3, 3}) == 24);
}

#include <random>
#include <map>
#include <limits>

TEST_CASE("DynamicNDimSegTree Stress Test") {
    auto min_func = [](int64_t a, int64_t b) { return std::min(a, b); };
    int64_t def = std::numeric_limits<int64_t>::max();
    vector<int> dims = {4, 4, 4, 4, 4};
    DynamicNDimSegTree<int64_t, decltype(min_func)> tree(dims, min_func, def);
    
    map<vector<int>, int64_t> grid_bf;
    for (int c0 = 0; c0 < 4; ++c0)
    for (int c1 = 0; c1 < 4; ++c1)
    for (int c2 = 0; c2 < 4; ++c2)
    for (int c3 = 0; c3 < 4; ++c3)
    for (int c4 = 0; c4 < 4; ++c4)
        grid_bf[{c0, c1, c2, c3, c4}] = def;

    mt19937 rng(42);
    uniform_int_distribution<int> dist_c(0, 3);
    uniform_int_distribution<int64_t> dist_v(1, 1000);

    for (int i = 0; i < 500; ++i) {
        vector<int> coords = {dist_c(rng), dist_c(rng), dist_c(rng), dist_c(rng), dist_c(rng)};
        int64_t val = dist_v(rng);
        tree.update(coords, val);
        grid_bf[coords] = val;
        
        vector<int> x_coords(5), y_coords(5);
        for (int d = 0; d < 5; ++d) {
            int a = dist_c(rng);
            int b = dist_c(rng);
            x_coords[d] = min(a, b);
            y_coords[d] = max(a, b);
        }

        int64_t expected = def;
        for (int c0 = x_coords[0]; c0 <= y_coords[0]; ++c0) {
            for (int c1 = x_coords[1]; c1 <= y_coords[1]; ++c1) {
                for (int c2 = x_coords[2]; c2 <= y_coords[2]; ++c2) {
                    for (int c3 = x_coords[3]; c3 <= y_coords[3]; ++c3) {
                        for (int c4 = x_coords[4]; c4 <= y_coords[4]; ++c4) {
                            expected = min(expected, grid_bf[{c0, c1, c2, c3, c4}]);
                        }
                    }
                }
            }
        }
        REQUIRE(tree.query_range(x_coords, y_coords) == expected);
    }
}
