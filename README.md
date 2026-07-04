# N-Dimensional Data Structures

A polyglot library providing generalized, highly optimized N-dimensional data structures for both Python and C++.

Whether you are working in Python for data science and backend engineering, or you need highly-optimized header-only C++ templates to copy-paste into competitive programming platforms like Codeforces and LeetCode, this repository has you covered.

## Data Structures Included

1. **Static N-Dim Prefix Sum** (`static_n_dim_prefix_sum`)
   Answers range sum queries offline. After all initial values are added, a `sweep()` computes prefix sums in $O(V)$ time, allowing subsequent $O(1)$ queries of any hyper-rectangular region.

2. **Static N-Dim Difference Array** (`static_n_dim_difference_array`)
   Applies many range additions offline. After deferred $O(1)$ updates are queued, a `sweep()` runs in $O(V)$ time to finalize the grid, allowing subsequent $O(1)$ point queries.

3. **Dynamic N-Dim Fenwick Tree** (`dynamic_n_dim_fenwick_tree`)
   Also known as a Binary Indexed Tree. Supports $O(\log^N(V))$ point updates and prefix sum queries across $N$ dimensions.

4. **Dynamic N-Dim Difference Fenwick Tree** (`dynamic_n_dim_diff_fenwick_tree`)
   A Fenwick tree modified to support $O(\log^N(V))$ range updates and point queries across $N$ dimensions.

5. **Dynamic N-Dim Range Fenwick Tree** (`dynamic_n_dim_range_fenwick_tree`)
   Supports $O(\log^N(V))$ range updates and range sum queries across $N$ dimensions. Internally uses $2^N$ algebraic Fenwick trees to correctly scale the range updates.

6. **Dynamic N-Dim Segment Tree** (`dynamic_n_dim_seg_tree`)
   Supports $O(\log^N(V))$ point updates and range queries across $N$ dimensions for any associative operation (e.g., `min`, `max`, `gcd`). Used when operations are not strictly invertible.

---

## Python Usage

The pure-Python package uses 1D list-flattening mathematics under the hood to ensure maximum performance and cache locality without relying on heavy C-extensions or Numpy.

### Installation
The package is officially available on PyPI:
```bash
pip install ndim_ds
```

### Example
```python
from ndim_ds import StaticNDimPrefixSum, DynamicNDimRangeFenwickTree

# --- Static Example ---
# Initialize a 10x10x10 static grid
grid = StaticNDimPrefixSum([10, 10, 10])

grid.add([1, 1, 1], 50)
grid.add([5, 5, 5], 25)
grid.sweep() # Finalize the grid in O(V) time

# Query the volume sum from (0, 0, 0) to (3, 3, 3) in O(1) time
total_static = grid.query_range([0, 0, 0], [3, 3, 3])

# --- Dynamic Example ---
# Initialize a 10x10x10 dynamic range tree
tree = DynamicNDimRangeFenwickTree([10, 10, 10])

# Add 50 to the bounding box from (1, 1, 1) to (5, 5, 5) in O(log^N(V)) time
tree.add_range([1, 1, 1], [5, 5, 5], 50)

# Query the volume sum from (0, 0, 0) to (3, 3, 3) in O(log^N(V)) time
total_dynamic = tree.query_range([0, 0, 0], [3, 3, 3])
```

---

## C++ (Competitive Programming) Usage

The `cpp/include/ndim/` directory contains highly optimized, header-only C++23 templates formatted for competitive programming. They use `#include <bits/stdc++.h>` and `using namespace std;`.

### Example
```cpp
#include "dynamic_n_dim_seg_tree.hpp"

// Example: 3D Segment Tree using std::min
auto min_func = [](int64_t a, int64_t b) { return min(a, b); };
int64_t def = 1e18;

// Initialize a 4x4x4 grid
DynamicNDimSegTree<int64_t, decltype(min_func)> tree({4, 4, 4}, min_func, def);

tree.update({1, 1, 1}, 5);
tree.update({2, 2, 2}, 10);

int64_t val = tree.query_range({0, 0, 0}, {3, 3, 3});
```

## Testing

Both ecosystems are rigorously tested with aggressive $10^3$ iteration stress tests across $5$-dimensional constraints.
- **Python**: Run `uv run pytest test/` from the `python/` directory.
- **C++**: Compile the files in `cpp/test/` using `g++ -std=c++23` and run the resulting executables.
