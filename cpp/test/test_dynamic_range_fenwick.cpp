#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"
#include "../include/ndim/dynamic_n_dim_range_fenwick_tree.hpp"
#include <vector>

using namespace ndim;

TEST_CASE("DynamicNDimRangeFenwickTree 2D Operations") {
    DynamicNDimRangeFenwickTree<int64_t> tree({5, 5});
    
    tree.add_range({1, 1}, {3, 3}, 10);
    tree.add_range({2, 2}, {4, 4}, 5);

    REQUIRE(tree.query_range({0, 0}, {0, 0}) == 0);
    
    REQUIRE(tree.query_range({1, 1}, {1, 1}) == 10);
    
    REQUIRE(tree.query_range({2, 2}, {3, 3}) == 60);
}

TEST_CASE("DynamicNDimRangeFenwickTree 3D Operations") {
    DynamicNDimRangeFenwickTree<int64_t> tree({4, 4, 4});
    
    tree.add_range({0, 0, 0}, {3, 3, 3}, 1);
    tree.add_range({1, 1, 1}, {2, 2, 2}, 2);

    REQUIRE(tree.query_range({0, 0, 0}, {0, 0, 0}) == 1);
    REQUIRE(tree.query_range({1, 1, 1}, {1, 1, 1}) == 3);
    
    REQUIRE(tree.query_range({1, 1, 1}, {2, 2, 2}) == 24);
    
    REQUIRE(tree.query_range({0, 0, 0}, {3, 3, 3}) == 80);
}

#include <random>
#include <map>

TEST_CASE("DynamicNDimRangeFenwickTree Stress Test") {
    vector<int> dims = {4, 4, 4, 4, 4};
    DynamicNDimRangeFenwickTree<int64_t> tree(dims);
    map<vector<int>, int64_t> grid_bf;

    mt19937 rng(42);
    uniform_int_distribution<int> dist_c(0, 3);
    uniform_int_distribution<int64_t> dist_v(1, 10);

    for (int i = 0; i < 500; ++i) {
        vector<int> x_coords(5), y_coords(5);
        for (int d = 0; d < 5; ++d) {
            int a = dist_c(rng);
            int b = dist_c(rng);
            x_coords[d] = min(a, b);
            y_coords[d] = max(a, b);
        }
        int64_t val = dist_v(rng);
        tree.add_range(x_coords, y_coords, val);
        
        for (int c0 = x_coords[0]; c0 <= y_coords[0]; ++c0) {
            for (int c1 = x_coords[1]; c1 <= y_coords[1]; ++c1) {
                for (int c2 = x_coords[2]; c2 <= y_coords[2]; ++c2) {
                    for (int c3 = x_coords[3]; c3 <= y_coords[3]; ++c3) {
                        for (int c4 = x_coords[4]; c4 <= y_coords[4]; ++c4) {
                            grid_bf[{c0, c1, c2, c3, c4}] += val;
                        }
                    }
                }
            }
        }
        
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
        REQUIRE(tree.query_range(x_coords, y_coords) == expected);
    }
}
