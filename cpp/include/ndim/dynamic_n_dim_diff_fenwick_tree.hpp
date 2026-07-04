#pragma once
#include "dynamic_n_dim_fenwick_tree.hpp"

namespace ndim {

template <typename T = long long>
class DynamicNDimDiffFenwickTree {
public:
    DynamicNDimDiffFenwickTree(const vector<int>& dims) : n_(dims.size()), tree_(dims), dims_(dims) {}

    void add_range(const vector<int>& x_coords, const vector<int>& y_coords, T val) {
        if (x_coords.size() != n_ || y_coords.size() != n_) {
            throw invalid_argument("Coordinates must match dimensions");
        }
        size_t num_masks = 1ULL << n_;
        for (size_t mask = 0; mask < num_masks; ++mask) {
            vector<int> target_coords(n_);
            int sign = 1;
            bool valid = true;
            for (size_t d = 0; d < n_; ++d) {
                if ((mask >> d) & 1) {
                    int c = y_coords[d] + 1;
                    if (c >= dims_[d]) {
                        valid = false;
                        break;
                    }
                    target_coords[d] = c;
                    sign = -sign;
                } else {
                    target_coords[d] = max(0, x_coords[d]);
                }
            }
            if (valid) {
                if (sign == 1) tree_.add(target_coords, val);
                else tree_.add(target_coords, -val);
            }
        }
    }

    T query_point(const vector<int>& coords) const {
        return tree_.query_prefix(coords);
    }

private:
    size_t n_;
    DynamicNDimFenwickTree<T> tree_;
    vector<int> dims_;
};

} // namespace ndim
