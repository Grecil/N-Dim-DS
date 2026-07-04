# N-Dim Data Structures

A collection of optimized N-dimensional data structures in Python, focused on memory efficiency and fast queries.

## Data Structures

The project contains several N-dimensional data structures optimized for different scenarios. They are located in the `src/` directory.

### Static Structures
- **Static N-Dim Prefix Sum (`static_n_dim_prefix_sum.py`)**: An N-dimensional grid optimized for offline point updates followed by rapid range queries.
- **Static N-Dim Difference Array (`static_n_dim_difference_array.py`)**: An N-dimensional grid initialized for bounding box range updates and later swept for querying specific points or regions.

### Dynamic Structures
- **Dynamic N-Dim Fenwick Tree (`dynamic_n_dim_fenwick_tree.py`)**: A Binary Indexed Tree supporting dynamic N-dimensional point updates and prefix/range queries in logarithmic time.
- **Dynamic N-Dim Diff Fenwick Tree (`dynamic_n_dim_diff_fenwick_tree.py`)**: A difference array over a Fenwick Tree, optimizing for dynamic N-dimensional range updates and exact point queries.
- **Dynamic N-Dim Range Fenwick Tree (`dynamic_n_dim_range_fenwick_tree.py`)**: A Fenwick Tree utilizing the 2^N algebraic expansion trick to support both dynamic range updates and range queries in logarithmic time.
- **Dynamic N-Dim Segment Tree (`dynamic_n_dim_seg_tree.py`)**: An iterative Segment Tree perfect for dynamic idempotent operations (Min, Max, GCD) across N-dimensional ranges.

## Directory Structure
- `src/`: Contains all the Python implementations of the data structures.
- `test/`: Intended directory for unit tests and verification logic.
