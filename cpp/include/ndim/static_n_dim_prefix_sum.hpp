#pragma once
#include <bits/stdc++.h>
using namespace std;

namespace ndim {

/**
 * A static N-dimensional Prefix Sum Array.
 * Supports O(1) deferred point updates across N dimensions.
 * Purpose: Efficiently answers range sum queries offline. After all initial values are added, a sweep() method
 * runs in O(V) time to compute prefix sums, allowing subsequent O(1) queries of any hyper-rectangular region.
 */
template <typename T = long long>
class StaticNDimPrefixSum {
public:
    StaticNDimPrefixSum(const vector<int>& dims) : dims_(dims), n_(dims.size()), is_swept_(false) {
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
     * @brief Sets or adds an initial value to a point.
     * @param coords 0-based point coordinates.
     * @param val Value to add to the point.
     */
    void add(const vector<int>& coords, T val) {
        if (is_swept_) {
            throw runtime_error("Cannot add points after the grid has been swept.");
        }
        if (coords.size() != n_) {
            throw invalid_argument("Coordinates must match the number of dimensions");
        }
        size_t idx = 0;
        for (size_t d = 0; d < n_; ++d) {
            if (coords[d] < 0 || coords[d] >= dims_[d]) {
                throw out_of_range("Coordinate out of bounds");
            }
            idx += coords[d] * strides_[d];
        }
        arr_[idx] += val;
    }

    /**
     * @brief Computes the prefix sums across all dimensions.
     * Must be called exactly once after all initial values are added, and before any queries.
     */
    void sweep() {
        for (size_t d = 0; d < n_; ++d) {
            size_t stride = strides_[d];
            size_t dim_size = dims_[d];
            for (size_t i = 0; i < total_size_; ++i) {
                if ((i / stride) % dim_size > 0) {
                    arr_[i] += arr_[i - stride];
                }
            }
        }
        is_swept_ = true;
    }

    /**
     * @brief Queries the total sum within a hyper-rectangular bounding box after sweeping.
     * @param x_coords 0-based starting indices (lower bounds).
     * @param y_coords 0-based ending indices (inclusive upper bounds).
     * @return T The sum within the range.
     */
    T query_range(const vector<int>& x_coords, const vector<int>& y_coords) const {
        if (!is_swept_) {
            throw runtime_error("Must call sweep() before querying.");
        }
        if (x_coords.size() != n_ || y_coords.size() != n_) {
            throw invalid_argument("Coordinates must match the number of dimensions");
        }
        T ans = 0;
        size_t num_masks = 1ULL << n_;
        for (size_t mask = 0; mask < num_masks; ++mask) {
            size_t idx = 0;
            int sign = 1;
            bool valid = true;
            for (size_t d = 0; d < n_; ++d) {
                if ((mask >> d) & 1) {
                    int c = x_coords[d] - 1;
                    if (c < 0) {
                        valid = false;
                        break;
                    }
                    idx += c * strides_[d];
                    sign = -sign;
                } else {
                    int c = min(y_coords[d], dims_[d] - 1);
                    idx += c * strides_[d];
                }
            }
            if (valid) {
                if (sign == 1) ans += arr_[idx];
                else ans -= arr_[idx];
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
    bool is_swept_;
};

} // namespace ndim
