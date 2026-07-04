#pragma once
#include <bits/stdc++.h>
using namespace std;

namespace ndim {

/**
 * A dynamic N-dimensional Segment Tree.
 * Supports O(log^N(V)) point updates and range queries across N dimensions.
 * Purpose: Efficiently updates individual elements and queries associative operations (e.g. min, max, gcd) over a hyper-rectangular region.
 * Used when operations are not strictly invertible (unlike Fenwick trees which require subtraction).
 */
template <typename T = long long, typename Func = function<T(T, T)>>
class DynamicNDimSegTree {
public:
    DynamicNDimSegTree(const vector<int>& dims, Func func, T default_val)
        : dims_(dims), n_(dims.size()), func_(func), default_(default_val) {
        if (dims.empty()) {
            throw invalid_argument("Dimensions cannot be empty");
        }
        sizes_.resize(n_);
        for (size_t i = 0; i < n_; ++i) {
            sizes_[i] = 2 * dims_[i];
        }

        strides_.resize(n_, 1);
        for (int i = n_ - 2; i >= 0; --i) {
            strides_[i] = strides_[i + 1] * sizes_[i + 1];
        }

        total_size_ = strides_[0] * sizes_[0];
        arr_.assign(total_size_, default_);
    }

    /**
     * @brief Updates the value at a specific point.
     * @param coords 0-based point coordinates.
     * @param val New value to set at the point.
     */
    void update(const vector<int>& coords, T val) {
        if (coords.size() != n_) {
            throw invalid_argument("Coordinates must match dimensions");
        }
        vector<int> leaf_coords(n_);
        for (size_t d = 0; d < n_; ++d) {
            leaf_coords[d] = coords[d] + dims_[d];
        }

        vector<vector<int>> paths(n_);
        for (size_t d = 0; d < n_; ++d) {
            int p = leaf_coords[d];
            while (p >= 1) {
                paths[d].push_back(p);
                p >>= 1;
            }
        }

        vector<size_t> current_indices(n_, 0);
        while (true) {
            size_t idx = 0;
            vector<int> p_tuple(n_);
            for (size_t d = 0; d < n_; ++d) {
                p_tuple[d] = paths[d][current_indices[d]];
                idx += p_tuple[d] * strides_[d];
            }

            int split_dim = -1;
            for (size_t d = 0; d < n_; ++d) {
                if (p_tuple[d] != leaf_coords[d]) {
                    split_dim = d;
                    break;
                }
            }

            if (split_dim == -1) {
                arr_[idx] = val;
            } else {
                size_t left_idx = idx + p_tuple[split_dim] * strides_[split_dim];
                size_t right_idx = left_idx + strides_[split_dim];
                arr_[idx] = func_(arr_[left_idx], arr_[right_idx]);
            }

            size_t d = 0;
            for (; d < n_; ++d) {
                current_indices[d]++;
                if (current_indices[d] < paths[d].size()) break;
                current_indices[d] = 0;
            }
            if (d == n_) break;
        }
    }

    /**
     * @brief Queries the aggregate function result over a hyper-rectangular bounding box.
     * @param x_coords 0-based starting indices (lower bounds).
     * @param y_coords 0-based ending indices (inclusive upper bounds).
     * @return T The result of the associative function over the range.
     */
    T query_range(const vector<int>& x_coords, const vector<int>& y_coords) const {
        if (x_coords.size() != n_ || y_coords.size() != n_) {
            throw invalid_argument("Coordinates must match dimensions");
        }

        vector<vector<size_t>> dim_nodes(n_);
        for (size_t d = 0; d < n_; ++d) {
            if (x_coords[d] > y_coords[d]) return default_;

            int l = x_coords[d] + dims_[d];
            int r = y_coords[d] + dims_[d] + 1;

            while (l < r) {
                if (l % 2 == 1) {
                    dim_nodes[d].push_back(l * strides_[d]);
                    l++;
                }
                if (r % 2 == 1) {
                    r--;
                    dim_nodes[d].push_back(r * strides_[d]);
                }
                l >>= 1;
                r >>= 1;
            }
            if (dim_nodes[d].empty()) return default_;
        }

        T res = default_;
        vector<size_t> current_indices(n_, 0);
        while (true) {
            size_t idx = 0;
            for (size_t d = 0; d < n_; ++d) {
                idx += dim_nodes[d][current_indices[d]];
            }
            res = func_(res, arr_[idx]);

            size_t d = 0;
            for (; d < n_; ++d) {
                current_indices[d]++;
                if (current_indices[d] < dim_nodes[d].size()) break;
                current_indices[d] = 0;
            }
            if (d == n_) break;
        }

        return res;
    }

private:
    vector<int> dims_;
    size_t n_;
    vector<int> sizes_;
    vector<size_t> strides_;
    size_t total_size_;
    vector<T> arr_;
    Func func_;
    T default_;
};

} // namespace ndim
