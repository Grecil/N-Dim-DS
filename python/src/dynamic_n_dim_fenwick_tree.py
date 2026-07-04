import itertools

class DynamicNDimFenwickTree:
    def __init__(self, dims):
        """
        Initializes an N-dimensional dynamic grid (Fenwick Tree / Binary Indexed Tree)
        capable of processing point updates and range queries in logarithmic time.

        The grid uses 0-based indexing and allocates memory strictly for the requested
        dimensions, ensuring no wasted padding space.

        Args:
            dims (iterable of int): The strict maximum size of the grid in each dimension.
                                    For example, dims=[5, 5] allows coordinates from (0,0) to (4,4).
        """
        self.n = len(dims)
        self.dims = tuple(dims)

        self.strides = [1] * self.n
        for i in range(self.n - 2, -1, -1):
            self.strides[i] = self.strides[i + 1] * self.dims[i + 1]

        self.total_size = self.strides[0] * self.dims[0]
        self.arr = [0] * self.total_size

    def add(self, coords, val):
        """
        Adds a value to a specific single point in the N-dimensional grid.
        This operation updates the internal tree structure dynamically.

        Args:
            coords (iterable of int): 0-based coordinate indices of the point to update.
                                      Must be of length N.
            val (int or float): The value to add to the specified point.
        """
        dim_indices = []
        for d in range(self.n):
            idx_list = []
            i = coords[d] + 1
            while i <= self.dims[d]:
                idx_list.append((i - 1) * self.strides[d])
                i += i & (-i)
            dim_indices.append(idx_list)

        for offsets in itertools.product(*dim_indices):
            self.arr[sum(offsets)] += val

    def query_prefix(self, coords):
        """
        Calculates the dynamic sum of all elements from the origin (0, 0, ...)
        up to and including the given coordinates.

        Args:
            coords (iterable of int): 0-based ending indices (inclusive).
                                      Must be of length N.

        Returns:
            int or float: The prefix sum from the origin to the specified coordinates.
        """
        dim_indices = []
        for d in range(self.n):
            idx_list = []
            i = coords[d] + 1
            i = min(i, self.dims[d])
            while i > 0:
                idx_list.append((i - 1) * self.strides[d])
                i -= i & (-i)
            dim_indices.append(idx_list)

        res = 0
        for offsets in itertools.product(*dim_indices):
            res += self.arr[sum(offsets)]

        return res

    def query_range(self, x_coords, y_coords):
        """
        Calculates the dynamic sum of all elements within an N-dimensional bounding box.
        Uses the Inclusion-Exclusion principle over the $2^N$ vertices of the hyper-rectangle.

        Args:
            x_coords (iterable of int): 0-based starting indices (lower bounds) of the bounding box.
            y_coords (iterable of int): 0-based ending indices (inclusive upper bounds) of the bounding box.

        Returns:
            int or float: The sum of the elements within the specified bounding box.
        """
        ans = 0
        for mask in range(1 << self.n):
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
