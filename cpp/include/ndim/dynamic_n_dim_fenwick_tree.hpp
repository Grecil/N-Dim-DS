#pragma once
#include <bits/stdc++.h>
using namespace std;

namespace ndim {

/**
 * A dynamic N-dimensional Fenwick Tree (Binary Indexed Tree).
 * Supports O(log^N(V)) point updates and prefix sum queries across N dimensions.
 * Purpose: Efficiently updates individual elements and computes the sum of elements in an N-dimensional prefix.
 */
template <typename T = long long>
class DynamicNDimFenwickTree {
public:
    DynamicNDimFenwickTree(const vector<int>& dims) : dims_(dims), n_(dims.size()) {
        if (dims.empty()) {
            throw invalid_argument("Dimensions cannot be empty");
        }
        strides_.resize(n_, 1);
        for (int i = n_ - 2; i >= 0; --i) {
            strides_[i] = strides_[i + 1] * dims_[i + 1];
        }
        total_size_ = strides_[0] * dims_[0];
        arr_.assign(total_size_, T{0});
    }

    /**
     * @brief Adds a value to a specific point.
     * @param coords 0-based point coordinates.
     * @param val Value to add.
     */
    void add(const vector<int>& coords, T val) {
        if (coords.size() != n_) {
            throw invalid_argument("Coordinates must match dimensions");
        }
        vector<vector<size_t>> dim_indices(n_);
        for (size_t d = 0; d < n_; ++d) {
            int i = coords[d] + 1;
            while (i <= dims_[d]) {
                dim_indices[d].push_back((i - 1) * strides_[d]);
                i += i & (-i);
            }
        }

        vector<size_t> current_indices(n_, 0);
        while (true) {
            size_t idx = 0;
            for (size_t d = 0; d < n_; ++d) {
                idx += dim_indices[d][current_indices[d]];
            }
            arr_[idx] += val;

            size_t d = 0;
            for (; d < n_; ++d) {
                current_indices[d]++;
                if (current_indices[d] < dim_indices[d].size()) break;
                current_indices[d] = 0;
            }
            if (d == n_) break;
        }
    }

    /**
     * @brief Queries the prefix sum from origin to the given point.
     * @param coords 0-based ending point of the prefix (inclusive).
     * @return T The prefix sum.
     */
    T query_prefix(const vector<int>& coords) const {
        vector<vector<size_t>> dim_indices(n_);
        for (size_t d = 0; d < n_; ++d) {
            int i = min(coords[d] + 1, dims_[d]);
            while (i > 0) {
                dim_indices[d].push_back((i - 1) * strides_[d]);
                i -= i & (-i);
            }
            if (dim_indices[d].empty()) return 0;
        }

        T res = 0;
        vector<size_t> current_indices(n_, 0);
        while (true) {
            size_t idx = 0;
            for (size_t d = 0; d < n_; ++d) {
                idx += dim_indices[d][current_indices[d]];
            }
            res += arr_[idx];

            size_t d = 0;
            for (; d < n_; ++d) {
                current_indices[d]++;
                if (current_indices[d] < dim_indices[d].size()) break;
                current_indices[d] = 0;
            }
            if (d == n_) break;
        }
        return res;
    }

    /**
     * @brief Queries the sum within a hyper-rectangular bounding box.
     * @param x_coords 0-based starting indices (lower bounds) of the bounding box.
     * @param y_coords 0-based ending indices (inclusive upper bounds) of the bounding box.
     * @return T The sum within the range.
     */
    T query_range(const vector<int>& x_coords, const vector<int>& y_coords) const {
        if (x_coords.size() != n_ || y_coords.size() != n_) {
            throw invalid_argument("Coordinates must match dimensions");
        }
        T ans = 0;
        size_t num_masks = 1ULL << n_;
        for (size_t mask = 0; mask < num_masks; ++mask) {
            vector<int> q_coords(n_);
            int sign = 1;
            bool valid = true;
            for (size_t d = 0; d < n_; ++d) {
                if ((mask >> d) & 1) {
                    int c = x_coords[d] - 1;
                    if (c < 0) {
                        valid = false;
                        break;
                    }
                    q_coords[d] = c;
                    sign = -sign;
                } else {
                    q_coords[d] = y_coords[d];
                }
            }
            if (valid) {
                if (sign == 1) ans += query_prefix(q_coords);
                else ans -= query_prefix(q_coords);
            }
        }
        return ans;
    }

private:
    vector<int> dims_;
    size_t n_;
    vector<size_t> strides_;
    size_t total_size_;
    vector<T> arr_;
};

} // namespace ndim
