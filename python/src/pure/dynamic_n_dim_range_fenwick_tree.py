import itertools


class DynamicNDimRangeFenwickTree:
    def __init__(self, dims):
        """
        Initializes an N-dimensional Fenwick Tree capable of processing BOTH
        dynamic range updates and dynamic range queries in logarithmic time.

        This relies on the "2^N Algebraic Expansion Trick", maintaining 2^N
        virtual trees packed tightly into a single array to compute spatial volumes.

        Args:
            dims (iterable of int): The strict maximum size of the grid in each dimension.
                                    Memory footprint is O(2^N * D_0 * D_1 * ... * D_n).
        """
        self.n = len(dims)
        self.dims = tuple(dims)
        self.num_masks = 1 << self.n

        self.strides = [1] * self.n
        for i in range(self.n - 2, -1, -1):
            self.strides[i] = self.strides[i + 1] * self.dims[i + 1]

        self.total_size = self.strides[0] * self.dims[0]
        # We store 2^N algebraic states per coordinate
        self.arr = [[0] * self.num_masks for _ in range(self.total_size)]

    def _add_point(self, coords, val):
        """
        Internal method. Updates all 2^N algebraic states at a specific point.
        """
        # Precompute the 2^N algebraic terms to store for this marker
        mask_vals = [0] * self.num_masks
        for mask in range(self.num_masks):
            v = val
            for d in range(self.n):
                if (mask >> d) & 1:
                    # -C_d (using 1-based coordinates for the math)
                    v *= -(coords[d] + 1)
            mask_vals[mask] = v

        dim_indices = []
        for d in range(self.n):
            idx_list = []
            i = coords[d] + 1
            while i <= self.dims[d]:
                idx_list.append((i - 1) * self.strides[d])
                i += i & (-i)
            dim_indices.append(idx_list)

        for offsets in itertools.product(*dim_indices):
            idx = sum(offsets)
            for mask in range(self.num_masks):
                self.arr[idx][mask] += mask_vals[mask]

    def add_range(self, x_coords, y_coords, val):
        """
        Dynamically adds a value to all elements within an N-dimensional bounding box.

        Args:
            x_coords (iterable of int): 0-based starting indices (lower bounds).
            y_coords (iterable of int): 0-based ending indices (inclusive upper bounds).
            val (int or float): The value to add to the specified volume.
        """
        for mask in range(self.num_masks):
            p_coords = [0] * self.n
            sign = 1
            valid = True

            for d in range(self.n):
                if (mask >> d) & 1:
                    c = y_coords[d] + 1
                    if c >= self.dims[d]:
                        valid = False
                        break
                    p_coords[d] = c
                    sign = -sign
                else:
                    p_coords[d] = x_coords[d]

            if valid:
                self._add_point(p_coords, sign * val)

    def query_prefix(self, coords):
        """
        Calculates the spatial volume from the origin up to the given coordinates
        using the 2^N algebraic expansion.
        """
        # Clamp coordinates to ensure we don't calculate volume outside the grid bounds
        clamped_coords = [min(coords[d], self.dims[d] - 1) for d in range(self.n)]

        dim_indices = []
        for d in range(self.n):
            idx_list = []
            i = clamped_coords[d] + 1
            while i > 0:
                idx_list.append((i - 1) * self.strides[d])
                i -= i & (-i)
            dim_indices.append(idx_list)

        # Accumulate the prefix sums for all 2^N virtual trees
        tree_sums = [0] * self.num_masks
        for offsets in itertools.product(*dim_indices):
            idx = sum(offsets)
            for mask in range(self.num_masks):
                tree_sums[mask] += self.arr[idx][mask]

        # Reconstruct the true prefix sum using the algebraic multipliers
        res = 0
        for mask in range(self.num_masks):
            multiplier = 1
            for d in range(self.n):
                if not ((mask >> d) & 1):
                    # (X_d + 1) in 1-based math translates to (clamped + 2)
                    multiplier *= clamped_coords[d] + 2
            res += multiplier * tree_sums[mask]

        return res

    def query_range(self, x_coords, y_coords):
        """
        Calculates the dynamic sum of all elements within an N-dimensional bounding box.

        Args:
            x_coords (iterable of int): 0-based starting indices.
            y_coords (iterable of int): 0-based ending indices (inclusive).

        Returns:
            int or float: The dynamic sum within the specified volume.
        """
        ans = 0
        for mask in range(self.num_masks):
            q_coords = [0] * self.n
            sign = 1
            valid = True

            for d in range(self.n):
                if (mask >> d) & 1:
                    c = x_coords[d] - 1
                    if c < 0:
                        valid = False
                        break
                    q_coords[d] = c
                    sign = -sign
                else:
                    q_coords[d] = y_coords[d]

            if valid:
                ans += sign * self.query_prefix(q_coords)

        return ans
