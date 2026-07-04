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

    REQUIRE(tree.query_range({0, 0, 0}, {1, 1, 1}) == 6); // gcd(12, 18) = 6
    REQUIRE(tree.query_range({1, 1, 1}, {2, 2, 2}) == 6); // gcd(18, 24) = 6
    REQUIRE(tree.query_range({0, 0, 0}, {2, 2, 2}) == 6); // gcd(12, 18, 24) = 6
    REQUIRE(tree.query_range({2, 2, 2}, {3, 3, 3}) == 24);
}
