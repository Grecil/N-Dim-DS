#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"
#include "../include/ndim/dynamic_n_dim_diff_fenwick_tree.hpp"
#include <vector>

using namespace ndim;

TEST_CASE("DynamicNDimDiffFenwickTree 2D Operations") {
    DynamicNDimDiffFenwickTree<int64_t> tree({5, 5});
    
    tree.add_range({1, 1}, {3, 3}, 10);
    tree.add_range({2, 2}, {4, 4}, 5);

    REQUIRE(tree.query_point({0, 0}) == 0);
    REQUIRE(tree.query_point({1, 1}) == 10);
    REQUIRE(tree.query_point({2, 2}) == 15);
    REQUIRE(tree.query_point({3, 3}) == 15);
    REQUIRE(tree.query_point({4, 4}) == 5);
}

TEST_CASE("DynamicNDimDiffFenwickTree 3D Operations") {
    DynamicNDimDiffFenwickTree<int64_t> tree({3, 3, 3});
    
    tree.add_range({0, 0, 0}, {2, 2, 2}, 1);
    tree.add_range({1, 1, 1}, {2, 2, 2}, 2);

    REQUIRE(tree.query_point({0, 0, 0}) == 1);
    REQUIRE(tree.query_point({1, 1, 1}) == 3);
    REQUIRE(tree.query_point({2, 2, 2}) == 3);
}
