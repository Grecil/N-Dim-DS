#pragma once
#include <bits/stdc++.h>
using namespace std;

namespace ndim {

template <typename T = int64_t>
class StaticNDimDifferenceArray {
public:
    StaticNDimDifferenceArray(const vector<int>& dims) : dims_(dims), n_(dims.size()), is_swept_(false) {
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

    void add_range(const vector<int>& x_coords, const vector<int>& y_coords, T val) {
        if (is_swept_) {
            throw runtime_error("Cannot add points after the grid has been swept.");
        }
        if (x_coords.size() != n_ || y_coords.size() != n_) {
            throw invalid_argument("Coordinates must match the number of dimensions");
        }
        size_t num_masks = 1ULL << n_;
        for (size_t mask = 0; mask < num_masks; ++mask) {
            size_t idx = 0;
            int sign = 1;
            bool valid = true;
            for (size_t d = 0; d < n_; ++d) {
                if ((mask >> d) & 1) {
                    int c = y_coords[d] + 1;
                    if (c >= dims_[d]) {
                        valid = false;
                        break;
                    }
                    idx += c * strides_[d];
                    sign = -sign;
                } else {
                    int c = max(0, x_coords[d]);
                    idx += c * strides_[d];
                }
            }
            if (valid) {
                if (sign == 1) arr_[idx] += val;
                else arr_[idx] -= val;
            }
        }
    }

    void sweep() {
        for (int iter = 0; iter < 2; ++iter) {
            for (size_t d = 0; d < n_; ++d) {
                size_t stride = strides_[d];
                size_t dim_size = dims_[d];
                for (size_t i = 0; i < total_size_; ++i) {
                    if ((i / stride) % dim_size > 0) {
                        arr_[i] += arr_[i - stride];
                    }
                }
            }
        }
        is_swept_ = true;
    }

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
