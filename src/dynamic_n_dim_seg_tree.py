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

        # Iterative seg trees require exactly 2 * D memory per dimension
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
        # Map 0-based coordinates to their exact leaf nodes in the iterative tree
        leaf_coords = [coords[d] + self.dims[d] for d in range(self.n)]

        # Generate the path from the leaf to the root for every dimension
        paths = []
        for d in range(self.n):
            path = []
            p = leaf_coords[d]
            while p >= 1:
                path.append(p)
                p >>= 1
            paths.append(path)

        # itertools.product inherently generates the nodes in bottom-up topological order
        for p_tuple in itertools.product(*paths):
            idx = 0
            for d in range(self.n):
                idx += p_tuple[d] * self.strides[d]

            # Find the first dimension where this node is a parent (not a leaf)
            split_dim = -1
            for d in range(self.n):
                if p_tuple[d] != leaf_coords[d]:
                    split_dim = d
                    break

            # If it's the exact N-dimensional leaf, set the value directly
            if split_dim == -1:
                self.arr[idx] = val
            # Otherwise, merge its two children along the split dimension
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

        # Step 1: For each dimension, find the exact set of 1D segment tree nodes
        # that perfectly cover the requested range [x, y].
        for d in range(self.n):
            if x_coords[d] > y_coords[d]:
                return self.default

            l = x_coords[d] + self.dims[d]
            r = y_coords[d] + self.dims[d] + 1
            nodes = []

            # Standard iterative 1D segment tree traversal
            while l < r:
                if l % 2 == 1:
                    nodes.append(l * self.strides[d])
                    l += 1
                if r % 2 == 1:
                    r -= 1
                    nodes.append(r * self.strides[d])
                l >>= 1
                r >>= 1

            dim_nodes.append(nodes)

        # Step 2: The exact N-dimensional nodes that cover our hyper-rectangle
        # are just the Cartesian product of the 1D nodes we found.
        res = self.default
        for offsets in itertools.product(*dim_nodes):
            res = self.func(res, self.arr[sum(offsets)])

        return res
