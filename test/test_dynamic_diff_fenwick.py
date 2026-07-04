import itertools
import random

from dynamic_n_dim_diff_fenwick_tree import DynamicNDimDiffFenwickTree


def test_dynamic_diff_fenwick_2d():
    tree = DynamicNDimDiffFenwickTree([5, 5])
    tree.add_range([1, 1], [2, 2], 5)
    assert tree.query_point([0, 0]) == 0
    assert tree.query_point([1, 1]) == 5
    assert tree.query_point([2, 2]) == 5
    assert tree.query_point([3, 3]) == 0
    assert tree.query_point([1, 2]) == 5
    tree.add_range([2, 2], [3, 3], 10)
    assert tree.query_point([1, 1]) == 5
    assert tree.query_point([2, 2]) == 15
    assert tree.query_point([3, 3]) == 10
    assert tree.query_point([4, 4]) == 0


def test_dynamic_diff_fenwick_3d():
    tree = DynamicNDimDiffFenwickTree([4, 4, 4])
    tree.add_range([0, 0, 0], [2, 2, 2], 100)
    assert tree.query_point([0, 0, 0]) == 100
    assert tree.query_point([2, 2, 2]) == 100
    assert tree.query_point([3, 3, 3]) == 0


def test_dynamic_diff_fenwick_stress():
    """Stress test with 5 dimensions (4x4x4x4x4 = 1024) and 1000 calls."""
    dims = [4, 4, 4, 4, 4]
    tree = DynamicNDimDiffFenwickTree(dims)
    grid_bf = {}

    random.seed(42)

    for _ in range(1000):
        x_coords = [random.randint(0, d - 1) for d in dims]
        y_coords = [random.randint(x, d - 1) for x, d in zip(x_coords, dims)]
        val = random.randint(1, 100)
        tree.add_range(x_coords, y_coords, val)

        ranges = [range(x, y + 1) for x, y in zip(x_coords, y_coords)]
        for c in itertools.product(*ranges):
            grid_bf[c] = grid_bf.get(c, 0) + val

    for _ in range(1000):
        coords = [random.randint(0, d - 1) for d in dims]
        expected = grid_bf.get(tuple(coords), 0)
        assert tree.query_point(coords) == expected
