#pragma once
#include <vector>
#include <stdexcept>
#include <cstdint>
#include <functional>
#include <algorithm>

namespace ndim {

template <typename T = int64_t, typename Func = std::function<T(T, T)>>
class DynamicNDimSegTree {
public:
    DynamicNDimSegTree(const std::vector<int>& dims, Func func, T default_val)
        : dims_(dims), n_(dims.size()), func_(func), default_(default_val) {
        if (dims.empty()) {
            throw std::invalid_argument("Dimensions cannot be empty");
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

    void update(const std::vector<int>& coords, T val) {
        if (coords.size() != n_) {
            throw std::invalid_argument("Coordinates must match dimensions");
        }
        std::vector<int> leaf_coords(n_);
        for (size_t d = 0; d < n_; ++d) {
            leaf_coords[d] = coords[d] + dims_[d];
        }

        std::vector<std::vector<int>> paths(n_);
        for (size_t d = 0; d < n_; ++d) {
            int p = leaf_coords[d];
            while (p >= 1) {
                paths[d].push_back(p);
                p >>= 1;
            }
        }

        std::vector<size_t> current_indices(n_, 0);
        while (true) {
            size_t idx = 0;
            std::vector<int> p_tuple(n_);
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

    T query_range(const std::vector<int>& x_coords, const std::vector<int>& y_coords) const {
        if (x_coords.size() != n_ || y_coords.size() != n_) {
            throw std::invalid_argument("Coordinates must match dimensions");
        }
        
        std::vector<std::vector<size_t>> dim_nodes(n_);
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
        std::vector<size_t> current_indices(n_, 0);
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
    std::vector<int> dims_;
    size_t n_;
    std::vector<int> sizes_;
    std::vector<size_t> strides_;
    size_t total_size_;
    std::vector<T> arr_;
    Func func_;
    T default_;
};

} // namespace ndim
