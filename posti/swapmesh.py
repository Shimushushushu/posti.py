import scipy.spatial


class SwapMesh:
    def __init__(self, x, U):
        self.x = x
        self.U = U
        self.dim = x.shape[1]
        self.tree = scipy.spatial.KDTree(x)
        self.max_d = 0

    def query(self, x):
        d, index = self.tree.query(x[0 : self.dim])
        # TODO(Shimushu): maybe linear interpolation
        self.max_d = max(d, self.max_d)
        return self.U[index]

    def max_distance(self):
        return self.max_d

    def reset(self):
        self.max_d = 0
