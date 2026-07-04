class StaticNDimPrefixSum:
    def __init__(self, dims):
        """
        Initializes an N-dimensional grid optimized for offline point updates 
        followed by rapid range queries using a static prefix sum array.
        
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
        self.is_swept = False

    def add(self, coords, val):
        """
        Adds a value to a specific point in the grid.
        This operation is deferred; the true prefix sum must be built by calling sweep()
        before any range queries can be made.
        
        Args:
            coords (iterable of int): 0-based coordinate indices of the point.
            val (int or float): The value to add to the specified point.
        """
        if self.is_swept:
            raise RuntimeError("Cannot add points after the grid has been swept.")
            
        idx = 0
        for d in range(self.n):
            idx += coords[d] * self.strides[d]
            
        self.arr[idx] += val

    def sweep(self):
        """
        Processes all accumulated point updates and computes the N-dimensional prefix sum.
        This transitions the data structure from an update phase to a query phase.
        Must be called exactly once after all additions and before any queries.
        """
        for d in range(self.n):
            stride = self.strides[d]
            dim_size = self.dims[d]
            for i in range(self.total_size):
                if (i // stride) % dim_size > 0:
                    self.arr[i] += self.arr[i - stride]
                    
        self.is_swept = True

    def query_range(self, x_coords, y_coords):
        """
        Calculates the exact sum of all elements within an N-dimensional bounding box.
        Uses the Inclusion-Exclusion principle over the $2^N$ vertices of the hyper-rectangle.
        
        Args:
            x_coords (iterable of int): 0-based starting indices (lower bounds) of the bounding box.
            y_coords (iterable of int): 0-based ending indices (inclusive upper bounds) of the bounding box.
            
        Returns:
            int or float: The sum of the elements within the specified bounding box.
        """
        if not self.is_swept:
            raise RuntimeError("Must call sweep() before querying.")
            
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
