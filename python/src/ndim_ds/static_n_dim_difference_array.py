class StaticNDimDifferenceArray:
    def __init__(self, dims):
        """
        Initializes an N-dimensional grid of zeros with the specified dimensions.

        Args:
            dims (iterable of int): The size of the grid in each dimension.
        """
        self.n = len(dims)
        self.dims = tuple(dims)

        self.strides = [1] * self.n
        for i in range(self.n - 2, -1, -1):
            self.strides[i] = self.strides[i + 1] * self.dims[i + 1]

        self.total_size = self.strides[0] * self.dims[0]
        self.arr = [0] * self.total_size
        self.is_swept = False

    def add_range(self, x_coords, y_coords, val):
        """
        Adds a value to all elements within an N-dimensional bounding box.
        This operation is deferred and will not be queryable until sweep() is called.

        Args:
            x_coords (iterable of int): 0-based starting indices of the bounding box
                for each dimension.
            y_coords (iterable of int): 0-based ending indices (inclusive) of the bounding box.
            val (int or float): The value to add to the specified range.
        """
        if self.is_swept:
            raise RuntimeError("Cannot add after sweep.")

        for mask in range(1 << self.n):
            idx = 0
            sign = 1
            valid = True
            for d in range(self.n):
                if (mask >> d) & 1:
                    c = y_coords[d] + 1
                    if c >= self.dims[d]:
                        valid = False
                        break
                    idx += c * self.strides[d]
                    sign = -sign
                else:
                    c = x_coords[d]
                    idx += c * self.strides[d]

            if valid:
                self.arr[idx] += sign * val

    def sweep(self):
        """
        Processes all accumulated range updates and prepares the grid for range queries.
        This transitions the data structure from an update phase to a query phase.
        Must be called exactly once after all additions and before any queries.
        """
        for _ in range(2):
            for d in range(self.n):
                stride = self.strides[d]
                dim_size = self.dims[d]
                for i in range(self.total_size):
                    if (i // stride) % dim_size > 0:
                        self.arr[i] += self.arr[i - stride]
        self.is_swept = True

    def query_range(self, x_coords, y_coords):
        """
        Calculates the sum of all elements within an N-dimensional bounding box.

        Args:
            x_coords (iterable of int): 0-based starting indices of the bounding box.
            y_coords (iterable of int): 0-based ending indices (inclusive) of the bounding box.

        Returns:
            int or float: The sum of the elements in the specified range.
        """
        if not self.is_swept:
            raise RuntimeError("Must sweep before query.")

        ans = 0
        for mask in range(1 << self.n):
            idx = 0
            sign = 1
            valid = True
            for d in range(self.n):
                if (mask >> d) & 1:
                    c = x_coords[d] - 1
                    if c < 0:
                        valid = False
                        break
                    idx += c * self.strides[d]
                    sign = -sign
                else:
                    c = y_coords[d]
                    idx += c * self.strides[d]

            if valid:
                ans += sign * self.arr[idx]
        return ans
