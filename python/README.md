# N-Dimensional Data Structures

A library providing generalized N-dimensional data structures for C++ and Python.

## Data Structures Included

1. **Static N-Dim Prefix Sum** (`static_n_dim_prefix_sum`)
2. **Static N-Dim Difference Array** (`static_n_dim_difference_array`)
3. **Dynamic N-Dim Fenwick Tree** (`dynamic_n_dim_fenwick_tree`)
4. **Dynamic N-Dim Difference Fenwick Tree** (`dynamic_n_dim_diff_fenwick_tree`)
5. **Dynamic N-Dim Range Fenwick Tree** (`dynamic_n_dim_range_fenwick_tree`)
6. **Dynamic N-Dim Segment Tree** (`dynamic_n_dim_seg_tree`)

---

## C++ Usage

The `cpp/include/ndim/` directory contains header-only C++23 templates formatted for competitive programming. They use `#include <bits/stdc++.h>` and `using namespace std;`.

### Example
```cpp
#include "dynamic_n_dim_seg_tree.hpp"

// Example: 3D Segment Tree using std::min
auto min_func = [](int64_t a, int64_t b) { return min(a, b); };
int64_t def = 1e18;

// 4x4x4 grid
DynamicNDimSegTree<int64_t, decltype(min_func)> tree({4, 4, 4}, min_func, def);

tree.update({1, 1, 1}, 5);
tree.update({2, 2, 2}, 10);

int64_t val = tree.query_range({0, 0, 0}, {3, 3, 3});
```

---

## Python Usage

The `python/src/` directory contains the Python implementations. The package is available on PyPI.

### Installation
```bash
pip install ndim_ds
```

### Example
```python
from static_n_dim_prefix_sum import StaticNDimPrefixSum

# 10x10x10 grid
grid = StaticNDimPrefixSum([10, 10, 10])

grid.add([1, 1, 1], 50)
grid.add([5, 5, 5], 25)
grid.sweep()

total = grid.query_range([0, 0, 0], [3, 3, 3])
```

## Testing

- **Python**: Run `uv run pytest test/`
- **C++**: Compile the files in `cpp/test/` using `g++ -std=c++23` and run the resulting executables.
