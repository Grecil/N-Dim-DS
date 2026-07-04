import itertools

class DynamicNDimDiffFenwickTree:
    def __init__(self, dims):
        """
        Initializes an N-dimensional dynamic difference array powered by a Fenwick Tree.
        This structure is optimized for applying values across continuous ranges
        and querying the exact value of any single point dynamically.

        The grid uses 0-based indexing and allocates memory strictly for the requested
        dimensions to prevent exponential memory bloat.

        Args:
            dims (iterable of int): The strict maximum size of the grid in each dimension.
                                    For example, dims=[10, 10, 10] creates a 3D grid
                                    allowing coordinates from (0,0,0) to (9,9,9).
        """
        self.n = len(dims)
        self.dims = tuple(dims)

        self.strides = [1] * self.n
        for i in range(self.n - 2, -1, -1):
            self.strides[i] = self.strides[i + 1] * self.dims[i + 1]

        self.total_size = self.strides[0] * self.dims[0]
        self.arr = [0] * self.total_size

    def _add_point(self, coords, val):
        """
        Internal helper method to apply a value to a specific point in the Fenwick Tree.
        Uses 1-based bitwise traversal mapped safely to the 0-based flat array.
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

    def add_range(self, x_coords, y_coords, val):
        """
        Dynamically adds a value to all elements within an N-dimensional bounding box.
        Places 2^N markers at the boundaries using the Inclusion-Exclusion principle.

        Args:
            x_coords (iterable of int): 0-based starting indices (lower bounds) of the bounding box.
            y_coords (iterable of int): 0-based ending indices (inclusive upper bounds).
            val (int or float): The value to add to every point inside the specified range.
        """
        for mask in range(1 << self.n):
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

    def query_point(self, coords):
        """
        Retrieves the exact, current value at a specific point in the N-dimensional grid.
        Since this tree stores a difference array, a point's value is calculated
        by summing all the difference markers from the origin up to the given coordinates.

        Args:
            coords (iterable of int): 0-based coordinates of the exact point to query.

        Returns:
            int or float: The current cumulative value at the specified point.
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
