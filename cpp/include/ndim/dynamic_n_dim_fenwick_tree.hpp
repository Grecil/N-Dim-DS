#pragma once
#include <vector>
#include <stdexcept>
#include <cstdint>
#include <algorithm>

namespace ndim {

template <typename T = int64_t>
class DynamicNDimFenwickTree {
public:
    DynamicNDimFenwickTree(const std::vector<int>& dims) : dims_(dims), n_(dims.size()) {
        if (dims.empty()) {
            throw std::invalid_argument("Dimensions cannot be empty");
        }
        strides_.resize(n_, 1);
        for (int i = n_ - 2; i >= 0; --i) {
            strides_[i] = strides_[i + 1] * dims_[i + 1];
        }
        total_size_ = strides_[0] * dims_[0];
        arr_.assign(total_size_, T{0});
    }

    void add(const std::vector<int>& coords, T val) {
        if (coords.size() != n_) {
            throw std::invalid_argument("Coordinates must match dimensions");
        }
        std::vector<std::vector<size_t>> dim_indices(n_);
        for (size_t d = 0; d < n_; ++d) {
            int i = coords[d] + 1;
            while (i <= dims_[d]) {
                dim_indices[d].push_back((i - 1) * strides_[d]);
                i += i & (-i);
            }
        }
        
        std::vector<size_t> current_indices(n_, 0);
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

    T query_prefix(const std::vector<int>& coords) const {
        std::vector<std::vector<size_t>> dim_indices(n_);
        for (size_t d = 0; d < n_; ++d) {
            int i = std::min(coords[d] + 1, dims_[d]);
            while (i > 0) {
                dim_indices[d].push_back((i - 1) * strides_[d]);
                i -= i & (-i);
            }
            if (dim_indices[d].empty()) return 0;
        }
        
        T res = 0;
        std::vector<size_t> current_indices(n_, 0);
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

    T query_range(const std::vector<int>& x_coords, const std::vector<int>& y_coords) const {
        if (x_coords.size() != n_ || y_coords.size() != n_) {
            throw std::invalid_argument("Coordinates must match dimensions");
        }
        T ans = 0;
        size_t num_masks = 1ULL << n_;
        for (size_t mask = 0; mask < num_masks; ++mask) {
            std::vector<int> q_coords(n_);
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
    std::vector<int> dims_;
    size_t n_;
    std::vector<size_t> strides_;
    size_t total_size_;
    std::vector<T> arr_;
};

} // namespace ndim
