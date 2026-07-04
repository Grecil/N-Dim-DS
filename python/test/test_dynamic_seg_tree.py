import itertools
import math
import random

from dynamic_n_dim_seg_tree import DynamicNDimSegTree


def test_dynamic_seg_tree_2d():
    tree = DynamicNDimSegTree([4, 4], func=min, default=float("inf"))
    tree.update([0, 0], 10)
    tree.update([1, 1], 5)
    tree.update([2, 2], 15)
    tree.update([3, 3], 20)
    assert tree.query_range([0, 0], [3, 3]) == 5
    assert tree.query_range([0, 0], [0, 0]) == 10
    assert tree.query_range([2, 2], [3, 3]) == 15
    assert tree.query_range([0, 2], [1, 3]) == float("inf")


def test_dynamic_seg_tree_gcd():
    tree = DynamicNDimSegTree([2, 2], func=math.gcd, default=0)
    tree.update([0, 0], 24)
    tree.update([1, 1], 36)
    assert tree.query_range([0, 0], [1, 1]) == 12
    assert tree.query_range([0, 0], [0, 0]) == 24


def test_dynamic_seg_tree_max():
    tree = DynamicNDimSegTree([3, 3, 3], func=max, default=float("-inf"))
    tree.update([0, 0, 0], 100)
    tree.update([1, 1, 1], 500)
    tree.update([2, 2, 2], 200)
    assert tree.query_range([0, 0, 0], [2, 2, 2]) == 500
    assert tree.query_range([0, 0, 0], [0, 1, 1]) == 100


def test_dynamic_seg_tree_stress():
    """Stress test with 5 dimensions (4x4x4x4x4 = 1024) and 1000 calls."""
    dims = [4, 4, 4, 4, 4]
    tree = DynamicNDimSegTree(dims, func=min, default=float("inf"))
    grid_bf = {}

    random.seed(42)

    for _ in range(1000):
        coords = [random.randint(0, d - 1) for d in dims]
        val = random.randint(1, 1000)
        tree.update(coords, val)
        c_tuple = tuple(coords)
        grid_bf[c_tuple] = val

    for _ in range(1000):
        x_coords = [random.randint(0, d - 1) for d in dims]
        y_coords = [random.randint(x, d - 1) for x, d in zip(x_coords, dims)]

        expected = float("inf")
        ranges = [range(x, y + 1) for x, y in zip(x_coords, y_coords)]
        for c in itertools.product(*ranges):
            expected = min(expected, grid_bf.get(c, float("inf")))

        assert tree.query_range(x_coords, y_coords) == expected
