#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"
#include "../include/ndim/dynamic_n_dim_fenwick_tree.hpp"
#include <vector>

using namespace ndim;

TEST_CASE("DynamicNDimFenwickTree 2D Operations") {
    DynamicNDimFenwickTree<int64_t> tree({5, 5});
    
    tree.add({1, 1}, 10);
    tree.add({2, 3}, 20);
    tree.add({4, 4}, 50);

    REQUIRE(tree.query_range({0, 0}, {0, 0}) == 0);
    REQUIRE(tree.query_range({1, 1}, {1, 1}) == 10);
    REQUIRE(tree.query_range({0, 0}, {2, 3}) == 30);
    REQUIRE(tree.query_range({2, 3}, {4, 4}) == 70);
}

TEST_CASE("DynamicNDimFenwickTree 3D Operations") {
    DynamicNDimFenwickTree<int64_t> tree({3, 3, 3});
    
    tree.add({0, 0, 0}, 1);
    tree.add({1, 1, 1}, 2);
    tree.add({2, 2, 2}, 3);

    REQUIRE(tree.query_range({0, 0, 0}, {0, 0, 0}) == 1);
    REQUIRE(tree.query_range({0, 0, 0}, {1, 1, 1}) == 3);
    REQUIRE(tree.query_range({1, 1, 1}, {2, 2, 2}) == 5);
    REQUIRE(tree.query_range({2, 2, 2}, {2, 2, 2}) == 3);
}
