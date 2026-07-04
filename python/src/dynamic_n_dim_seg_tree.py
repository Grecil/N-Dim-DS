import itertools


class DynamicNDimSegTree:
    def __init__(self, dims, func=min, default=float("inf")):
        """
        Initializes an N-dimensional Iterative Segment Tree.
        Perfect for dynamic idempotent operations (Min, Max, GCD).

        The grid uses 0-based indexing. It leverages an iterative array layout
        to strictly bound memory to exactly 2^N * D_1 * D_2 ... * D_n.

        Args:
            dims (iterable of int): The strict maximum size of the grid in each dimension.
            func (callable): The function to merge two nodes (e.g., min, max, math.gcd).
            default (int or float): The identity value for the function.
        """
        self.n = len(dims)
        self.dims = tuple(dims)
        self.func = func
        self.default = default

        self.sizes = tuple(2 * d for d in dims)

        self.strides = [1] * self.n
        for i in range(self.n - 2, -1, -1):
            self.strides[i] = self.strides[i + 1] * self.sizes[i + 1]

        self.total_size = self.strides[0] * self.sizes[0]
        self.arr = [self.default] * self.total_size

    def update(self, coords, val):
        """
        Updates a specific point and iteratively recomputes the affected tree hierarchy.

        Args:
            coords (iterable of int): 0-based coordinate indices of the point.
            val (int or float): The new value for the specified point.
        """
        leaf_coords = [coords[d] + self.dims[d] for d in range(self.n)]

        paths = []
        for d in range(self.n):
            path = []
            p = leaf_coords[d]
            while p >= 1:
                path.append(p)
                p >>= 1
            paths.append(path)

        for p_tuple in itertools.product(*paths):
            idx = 0
            for d in range(self.n):
                idx += p_tuple[d] * self.strides[d]

            split_dim = -1
            for d in range(self.n):
                if p_tuple[d] != leaf_coords[d]:
                    split_dim = d
                    break

            if split_dim == -1:
                self.arr[idx] = val
            else:
                left_idx = idx + p_tuple[split_dim] * self.strides[split_dim]
                right_idx = left_idx + self.strides[split_dim]
                self.arr[idx] = self.func(self.arr[left_idx], self.arr[right_idx])

    def query_range(self, x_coords, y_coords):
        """
        Iteratively calculates the aggregate function over an N-dimensional bounding box.

        Args:
            x_coords (iterable of int): 0-based starting indices (lower bounds).
            y_coords (iterable of int): 0-based ending indices (inclusive upper bounds).

        Returns:
            int or float: The aggregated result (e.g., minimum value) in the volume.
        """
        dim_nodes = []

        for d in range(self.n):
            if x_coords[d] > y_coords[d]:
                return self.default

            left_idx = x_coords[d] + self.dims[d]
            r = y_coords[d] + self.dims[d] + 1
            nodes = []

            while left_idx < r:
                if left_idx % 2 == 1:
                    nodes.append(left_idx * self.strides[d])
                    left_idx += 1
                if r % 2 == 1:
                    r -= 1
                    nodes.append(r * self.strides[d])
                left_idx >>= 1
                r >>= 1

            dim_nodes.append(nodes)

        res = self.default
        for offsets in itertools.product(*dim_nodes):
            res = self.func(res, self.arr[sum(offsets)])

        return res
