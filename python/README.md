# N-Dimensional Data Structures (Python)

A pure-Python library providing generalized, highly optimized N-dimensional data structures. These structures use 1D list-flattening mathematics under the hood to ensure maximum performance and cache locality without relying on heavy C-extensions or Numpy.

## Installation

The package is officially available on PyPI:

```bash
pip install ndim_ds
```

## Data Structures Included

1. **Static N-Dim Prefix Sum** (`StaticNDimPrefixSum`)
   Answers range sum queries offline. After all initial values are added, a `sweep()` computes prefix sums in $O(V)$ time, allowing subsequent $O(1)$ queries of any hyper-rectangular region.

2. **Static N-Dim Difference Array** (`StaticNDimDifferenceArray`)
   Applies many range additions offline. After deferred $O(1)$ updates are queued, a `sweep()` runs in $O(V)$ time to finalize the grid, allowing subsequent $O(1)$ point queries.

3. **Dynamic N-Dim Fenwick Tree** (`DynamicNDimFenwickTree`)
   Also known as a Binary Indexed Tree. Supports $O(\log^N(V))$ point updates and prefix sum queries across $N$ dimensions.

4. **Dynamic N-Dim Difference Fenwick Tree** (`DynamicNDimDiffFenwickTree`)
   A Fenwick tree modified to support $O(\log^N(V))$ range updates and point queries across $N$ dimensions.

5. **Dynamic N-Dim Range Fenwick Tree** (`DynamicNDimRangeFenwickTree`)
   Supports $O(\log^N(V))$ range updates and range sum queries across $N$ dimensions. Internally uses $2^N$ algebraic Fenwick trees to correctly scale the range updates.

6. **Dynamic N-Dim Segment Tree** (`DynamicNDimSegTree`)
   Supports $O(\log^N(V))$ point updates and range queries across $N$ dimensions for any associative operation (e.g., `min`, `max`, `gcd`).

---

## Global Usage

Once installed globally via `pip`, you can directly import the data structures into any Python script on your machine.

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

## Testing

The codebase is rigorously tested with aggressive $10^3$ iteration stress tests across $5$-dimensional constraints using `pytest`.
To run tests locally from source:
```bash
uv run pytest test/
```
