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
    
    // (1,1) -> 10
    REQUIRE(tree.query_range({1, 1}, {1, 1}) == 10);
    
    // (2,2) to (3,3) overlap both -> 15 each. (2x2 area = 4 points * 15 = 60)
    REQUIRE(tree.query_range({2, 2}, {3, 3}) == 60);
}

TEST_CASE("DynamicNDimRangeFenwickTree 3D Operations") {
    DynamicNDimRangeFenwickTree<int64_t> tree({4, 4, 4});
    
    tree.add_range({0, 0, 0}, {3, 3, 3}, 1);
    tree.add_range({1, 1, 1}, {2, 2, 2}, 2);

    REQUIRE(tree.query_range({0, 0, 0}, {0, 0, 0}) == 1);
    REQUIRE(tree.query_range({1, 1, 1}, {1, 1, 1}) == 3);
    
    // (1,1,1) to (2,2,2) is a 2x2x2 block (8 points), each has value 3. 8*3 = 24
    REQUIRE(tree.query_range({1, 1, 1}, {2, 2, 2}) == 24);
    
    // whole grid is 64 points * 1 = 64 + inner 8 points * 2 = 16 => 80
    REQUIRE(tree.query_range({0, 0, 0}, {3, 3, 3}) == 80);
}
