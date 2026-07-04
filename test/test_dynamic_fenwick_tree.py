import itertools
import random

from dynamic_n_dim_fenwick_tree import DynamicNDimFenwickTree


def test_dynamic_fenwick_tree_2d():
    tree = DynamicNDimFenwickTree([5, 5])
    tree.add([1, 1], 10)
    tree.add([2, 2], 5)
    tree.add([3, 3], 15)
    assert tree.query_range([1, 1], [1, 1]) == 10
    assert tree.query_range([0, 0], [2, 2]) == 15
    assert tree.query_range([1, 1], [4, 4]) == 30
    assert tree.query_range([2, 2], [3, 3]) == 20
    assert tree.query_range([0, 0], [0, 0]) == 0
    tree.add([2, 2], 5)
    assert tree.query_range([2, 2], [3, 3]) == 25


def test_dynamic_fenwick_tree_3d():
    tree = DynamicNDimFenwickTree([3, 3, 3])
    tree.add([0, 0, 0], 1)
    tree.add([1, 1, 1], 2)
    tree.add([2, 2, 2], 3)
    assert tree.query_range([0, 0, 0], [2, 2, 2]) == 6
    assert tree.query_range([0, 0, 0], [1, 1, 1]) == 3
    assert tree.query_range([1, 1, 1], [2, 2, 2]) == 5
    assert tree.query_prefix([2, 2, 2]) == 6


def test_dynamic_fenwick_tree_stress():
    """Stress test with 5 dimensions (4x4x4x4x4 = 1024) and 1000 calls."""
    dims = [4, 4, 4, 4, 4]
    tree = DynamicNDimFenwickTree(dims)
    grid_bf = {}

    random.seed(42)

    for _ in range(1000):
        coords = [random.randint(0, d - 1) for d in dims]
        val = random.randint(1, 100)
        tree.add(coords, val)
        c_tuple = tuple(coords)
        grid_bf[c_tuple] = grid_bf.get(c_tuple, 0) + val

    for _ in range(1000):
        x_coords = [random.randint(0, d - 1) for d in dims]
        y_coords = [random.randint(x, d - 1) for x, d in zip(x_coords, dims)]

        expected = 0
        ranges = [range(x, y + 1) for x, y in zip(x_coords, y_coords)]
        for c in itertools.product(*ranges):
            expected += grid_bf.get(c, 0)

        assert tree.query_range(x_coords, y_coords) == expected
