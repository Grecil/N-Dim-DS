import itertools
import random

import pytest

from pure.static_n_dim_difference_array import StaticNDimDifferenceArray


def test_static_difference_array_2d():
    grid = StaticNDimDifferenceArray([4, 4])
    grid.add_range([1, 1], [2, 2], 5)
    grid.add_range([0, 0], [1, 2], 3)
    grid.sweep()
    assert grid.query_range([0, 0], [0, 0]) == 3
    assert grid.query_range([1, 1], [1, 1]) == 8
    assert grid.query_range([2, 2], [2, 2]) == 5
    assert grid.query_range([3, 3], [3, 3]) == 0
    assert grid.query_range([1, 1], [2, 2]) == 26


def test_static_difference_array_3d():
    grid = StaticNDimDifferenceArray([3, 3, 3])
    grid.add_range([0, 0, 0], [1, 1, 1], 10)
    grid.sweep()
    assert grid.query_range([0, 0, 0], [0, 0, 0]) == 10
    assert grid.query_range([1, 1, 1], [1, 1, 1]) == 10
    assert grid.query_range([2, 2, 2], [2, 2, 2]) == 0
    assert grid.query_range([0, 0, 0], [1, 1, 1]) == 80


def test_static_difference_array_errors():
    grid = StaticNDimDifferenceArray([2, 2])
    grid.add_range([0, 0], [1, 1], 1)
    with pytest.raises(RuntimeError):
        grid.query_range([0, 0], [1, 1])
    grid.sweep()
    with pytest.raises(RuntimeError):
        grid.add_range([0, 0], [0, 0], 5)


def test_static_difference_array_stress():
    """Stress test with 5 dimensions (4x4x4x4x4 = 1024) and 1000 calls."""
    dims = [4, 4, 4, 4, 4]
    grid = StaticNDimDifferenceArray(dims)
    grid_bf = {}

    random.seed(42)

    for _ in range(1000):
        x_coords = [random.randint(0, d - 1) for d in dims]
        y_coords = [random.randint(x, d - 1) for x, d in zip(x_coords, dims)]
        val = random.randint(1, 100)
        grid.add_range(x_coords, y_coords, val)

        ranges = [range(x, y + 1) for x, y in zip(x_coords, y_coords)]
        for c in itertools.product(*ranges):
            grid_bf[c] = grid_bf.get(c, 0) + val

    grid.sweep()

    for _ in range(1000):
        x_coords = [random.randint(0, d - 1) for d in dims]
        y_coords = [random.randint(x, d - 1) for x, d in zip(x_coords, dims)]

        expected = 0
        ranges = [range(x, y + 1) for x, y in zip(x_coords, y_coords)]
        for c in itertools.product(*ranges):
            expected += grid_bf.get(c, 0)

        assert grid.query_range(x_coords, y_coords) == expected
