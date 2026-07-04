# N-Dim Data Structures

A polyglot collection of optimized N-dimensional data structures, focused on memory efficiency and fast queries.

## Architecture & Ecosystems

This repository is designed as a **Polyglot Monorepo** to cater to three distinct development ecosystems:
1. **Pure Python** (General development and scripting)
2. **NumPy Python** (Data science and heavy matrix computations) - *Coming soon*
3. **C++ Templates** (Competitive programming and native performance) - *Coming soon*

## Directory Structure

- `python/`: Contains all Python implementations managed by `uv`.
  - `src/pure/`: Native Python implementations relying on standard libraries.
  - `src/numpy_impl/`: NumPy-accelerated implementations (planned).
  - `test/`: Comprehensive `pytest` suite including 5D stress tests.
- `cpp/`: Contains C++ implementations (planned).
  - `include/`: Header-only (`.hpp`) recursive template algorithms for easy copy-pasting.
  - `test/`: Catch2 unit tests via CMake.

## Data Structures (Python Pure)

The following algorithms are currently available in the `python/src/pure/` directory:

### Static Structures
- **Static N-Dim Prefix Sum**: An N-dimensional grid optimized for offline point updates followed by rapid range queries.
- **Static N-Dim Difference Array**: An N-dimensional grid initialized for bounding box range updates and later swept for querying specific points or regions.

### Dynamic Structures
- **Dynamic N-Dim Fenwick Tree**: A Binary Indexed Tree supporting dynamic N-dimensional point updates and prefix/range queries in logarithmic time.
- **Dynamic N-Dim Diff Fenwick Tree**: A difference array over a Fenwick Tree, optimizing for dynamic N-dimensional range updates and exact point queries.
- **Dynamic N-Dim Range Fenwick Tree**: A Fenwick Tree utilizing the $2^N$ algebraic expansion trick to support both dynamic range updates and range queries in logarithmic time.
- **Dynamic N-Dim Segment Tree**: An iterative Segment Tree perfect for dynamic idempotent operations (Min, Max, GCD) across N-dimensional ranges.

## Testing
To run the extensive stress test suite for the Python implementations:
```bash
cd python
uv run pytest test/ -v
```
