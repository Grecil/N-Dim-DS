# N-Dimensional Data Structures

A blazing-fast, polyglot (Python & C++) library providing generalized $N$-dimensional data structures. 

Whether you are working in Python for data science and backend engineering, or you need highly-optimized header-only C++ templates to copy-paste into competitive programming platforms like Codeforces and LeetCode, this repository has you covered.

## Data Structures Included
All structures support $O(1)$ or $O(\log^N(V))$ operations scaled effortlessly across any arbitrary number of dimensions.

1. **Static N-Dim Prefix Sum** (`static_n_dim_prefix_sum`)
2. **Static N-Dim Difference Array** (`static_n_dim_difference_array`)
3. **Dynamic N-Dim Fenwick Tree** (`dynamic_n_dim_fenwick_tree`)
4. **Dynamic N-Dim Difference Fenwick Tree** (`dynamic_n_dim_diff_fenwick_tree`)
5. **Dynamic N-Dim Range Fenwick Tree** (`dynamic_n_dim_range_fenwick_tree`)
6. **Dynamic N-Dim Segment Tree** (`dynamic_n_dim_seg_tree`)

---

## C++ (Competitive Programming)

The `cpp/include/ndim/` directory contains highly optimized, header-only C++17 templates designed specifically for **competitive programming**.

They use standard `<bits/stdc++.h>` format with `using namespace std;` to ensure they are instantly copy-pasteable into online judges without causing scope or header errors. Memory is maintained via flat, contiguous 1D vectors mapped algebraically to $N$ dimensions, preventing memory scattering and ensuring maximum CPU cache efficiency.

### Example Usage
```cpp
#include "dynamic_n_dim_seg_tree.hpp"

// Example: 3D Segment Tree using std::min
auto min_func = [](int64_t a, int64_t b) { return min(a, b); };
int64_t def = 1e18; // infinity

// 4x4x4 grid
DynamicNDimSegTree<int64_t, decltype(min_func)> tree({4, 4, 4}, min_func, def);

tree.update({1, 1, 1}, 5);
tree.update({2, 2, 2}, 10);

int64_t val = tree.query_range({0, 0, 0}, {3, 3, 3}); // 5
```

---

## Python (Data Science / General)

The `python/src/` directory contains the pure-Python implementations. They use 1D list-flattening mathematics similar to numpy under the hood, but operate purely on standard library primitives, ensuring maximum portability without heavy C-extension dependencies.

### Installation
The library is configured using `uv` via `pyproject.toml`.
```bash
cd python
uv sync
```

### Example Usage
```python
from dynamic_n_dim_range_fenwick_tree import DynamicNDimRangeFenwickTree

# 10x10x10 grid
tree = DynamicNDimRangeFenwickTree([10, 10, 10])

# Add 50 to the bounding box from (1, 1, 1) to (5, 5, 5)
tree.add_range([1, 1, 1], [5, 5, 5], 50)

# Query the volume sum from (0, 0, 0) to (3, 3, 3)
total = tree.query_range([0, 0, 0], [3, 3, 3])
```

## Testing

Both ecosystems are rigorously tested with aggressive $10^3$ iteration stress tests across $5$-dimensional constraints.
- **Python**: Uses `pytest`. Run `uv run pytest test/`
- **C++**: Uses `doctest`. Compile the files in `cpp/test/` using `g++ -std=c++17` and run the resulting executables.
