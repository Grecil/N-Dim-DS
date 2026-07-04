import itertools
import random

import pytest

from pure.static_n_dim_prefix_sum import StaticNDimPrefixSum


def test_static_prefix_sum_2d():
    grid = StaticNDimPrefixSum([3, 3])
    grid.add([0, 0], 1)
    grid.add([1, 1], 2)
    grid.add([2, 2], 3)
    grid.sweep()
    assert grid.query_range([0, 0], [0, 0]) == 1
    assert grid.query_range([1, 1], [1, 1]) == 2
    assert grid.query_range([2, 2], [2, 2]) == 3
    assert grid.query_range([0, 0], [2, 2]) == 6
    assert grid.query_range([0, 0], [1, 1]) == 3
    assert grid.query_range([0, 0], [0, 1]) == 1


def test_static_prefix_sum_3d():
    grid = StaticNDimPrefixSum([2, 2, 2])
    grid.add([0, 0, 0], 5)
    grid.add([1, 1, 1], 10)
    grid.sweep()
    assert grid.query_range([0, 0, 0], [1, 1, 1]) == 15
    assert grid.query_range([0, 0, 0], [0, 1, 1]) == 5
    assert grid.query_range([1, 1, 1], [1, 1, 1]) == 10


def test_static_prefix_sum_errors():
    grid = StaticNDimPrefixSum([2, 2])
    grid.add([0, 0], 1)
    with pytest.raises(RuntimeError):
        grid.query_range([0, 0], [1, 1])
    grid.sweep()
    with pytest.raises(RuntimeError):
        grid.add([1, 1], 5)


def test_static_prefix_sum_stress():
    """Stress test with 5 dimensions (4x4x4x4x4 = 1024) and 1000 calls."""
    dims = [4, 4, 4, 4, 4]
    grid = StaticNDimPrefixSum(dims)
    grid_bf = {}

    random.seed(42)

    for _ in range(1000):
        coords = [random.randint(0, d - 1) for d in dims]
        val = random.randint(1, 100)
        grid.add(coords, val)
        c_tuple = tuple(coords)
        grid_bf[c_tuple] = grid_bf.get(c_tuple, 0) + val

    grid.sweep()

    for _ in range(1000):
        x_coords = [random.randint(0, d - 1) for d in dims]
        y_coords = [random.randint(x, d - 1) for x, d in zip(x_coords, dims)]

        expected = 0
        ranges = [range(x, y + 1) for x, y in zip(x_coords, y_coords)]
        for c in itertools.product(*ranges):
            expected += grid_bf.get(c, 0)

        assert grid.query_range(x_coords, y_coords) == expected
