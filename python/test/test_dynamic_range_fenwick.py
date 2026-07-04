import itertools
import random

from pure.dynamic_n_dim_range_fenwick_tree import DynamicNDimRangeFenwickTree


def test_dynamic_range_fenwick_2d():
    tree = DynamicNDimRangeFenwickTree([4, 4])
    tree.add_range([0, 0], [1, 1], 5)
    assert tree.query_range([0, 0], [0, 0]) == 5
    assert tree.query_range([1, 1], [1, 1]) == 5
    assert tree.query_range([0, 0], [1, 1]) == 20
    assert tree.query_range([0, 0], [3, 3]) == 20
    tree.add_range([1, 1], [2, 2], 10)
    assert tree.query_range([1, 1], [1, 1]) == 15
    assert tree.query_range([0, 0], [3, 3]) == 60
    assert tree.query_prefix([1, 1]) == 30


def test_dynamic_range_fenwick_3d():
    tree = DynamicNDimRangeFenwickTree([3, 3, 3])
    tree.add_range([0, 0, 0], [1, 1, 1], 1)
    assert tree.query_range([0, 0, 0], [0, 0, 0]) == 1
    assert tree.query_range([0, 0, 0], [1, 1, 1]) == 8
    assert tree.query_range([0, 0, 0], [2, 2, 2]) == 8


def test_dynamic_range_fenwick_stress():
    """Stress test with 5 dimensions (4x4x4x4x4 = 1024) and 1000 calls."""
    dims = [4, 4, 4, 4, 4]
    tree = DynamicNDimRangeFenwickTree(dims)
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
        x_coords = [random.randint(0, d - 1) for d in dims]
        y_coords = [random.randint(x, d - 1) for x, d in zip(x_coords, dims)]

        expected = 0
        ranges = [range(x, y + 1) for x, y in zip(x_coords, y_coords)]
        for c in itertools.product(*ranges):
            expected += grid_bf.get(c, 0)

        assert tree.query_range(x_coords, y_coords) == expected
