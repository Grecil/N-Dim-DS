#pragma once
#include <bits/stdc++.h>
using namespace std;

namespace ndim {

template <typename T = long long>
class DynamicNDimRangeFenwickTree {
public:
    DynamicNDimRangeFenwickTree(const vector<int>& dims) : dims_(dims), n_(dims.size()) {
        if (dims.empty()) {
            throw invalid_argument("Dimensions cannot be empty");
        }
        strides_.resize(n_, 1);
        for (int i = n_ - 2; i >= 0; --i) {
            strides_[i] = strides_[i + 1] * dims_[i + 1];
        }
        total_size_ = strides_[0] * dims_[0];

        size_t num_trees = 1ULL << n_;
        trees_.assign(num_trees, vector<T>(total_size_, T{0}));
    }

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
                if (sign == 1) _point_add(target_coords, val);
                else _point_add(target_coords, -val);
            }
        }
    }

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
                if (sign == 1) ans += _prefix_query(q_coords);
                else ans -= _prefix_query(q_coords);
            }
        }
        return ans;
    }

private:
    vector<int> dims_;
    size_t n_;
    vector<size_t> strides_;
    size_t total_size_;
    vector<vector<T>> trees_;

    void _point_add(const vector<int>& coords, T val) {
        vector<vector<size_t>> dim_indices(n_);
        for (size_t d = 0; d < n_; ++d) {
            int i = coords[d] + 1;
            while (i <= dims_[d]) {
                dim_indices[d].push_back((i - 1) * strides_[d]);
                i += i & (-i);
            }
        }

        size_t num_masks = 1ULL << n_;
        vector<T> mask_vals(num_masks);
        for (size_t mask = 0; mask < num_masks; ++mask) {
            T v = val;
            for (size_t d = 0; d < n_; ++d) {
                if ((mask >> d) & 1) {
                    v *= -(coords[d] + 1);
                }
            }
            mask_vals[mask] = v;
        }

        for (size_t mask = 0; mask < num_masks; ++mask) {
            vector<size_t> current_indices(n_, 0);
            while (true) {
                size_t idx = 0;
                for (size_t d = 0; d < n_; ++d) {
                    idx += dim_indices[d][current_indices[d]];
                }
                trees_[mask][idx] += mask_vals[mask];

                size_t d = 0;
                for (; d < n_; ++d) {
                    current_indices[d]++;
                    if (current_indices[d] < dim_indices[d].size()) break;
                    current_indices[d] = 0;
                }
                if (d == n_) break;
            }
        }
    }

    T _prefix_query(const vector<int>& coords) const {
        vector<vector<size_t>> dim_indices(n_);
        for (size_t d = 0; d < n_; ++d) {
            int i = min(coords[d] + 1, dims_[d]);
            while (i > 0) {
                dim_indices[d].push_back((i - 1) * strides_[d]);
                i -= i & (-i);
            }
            if (dim_indices[d].empty()) return 0;
        }

        T total = 0;
        size_t num_masks = 1ULL << n_;
        for (size_t mask = 0; mask < num_masks; ++mask) {
            T sum_val = 0;
            vector<size_t> current_indices(n_, 0);
            while (true) {
                size_t idx = 0;
                for (size_t d = 0; d < n_; ++d) {
                    idx += dim_indices[d][current_indices[d]];
                }
                sum_val += trees_[mask][idx];

                size_t d = 0;
                for (; d < n_; ++d) {
                    current_indices[d]++;
                    if (current_indices[d] < dim_indices[d].size()) break;
                    current_indices[d] = 0;
                }
                if (d == n_) break;
            }

            T mult = 1;
            for (size_t d = 0; d < n_; ++d) {
                if (!((mask >> d) & 1)) {
                    mult *= (coords[d] + 2);
                }
            }
            total += sum_val * mult;
        }
        return total;
    }
};

} // namespace ndim
