import scipy.spatial
from scipy.interpolate import RBFInterpolator
import numpy as np


class SwapMesh:
    def __init__(self, x, U):
        self.x = x
        self.U = U
        self.dim = x.shape[1]
        self.tree = scipy.spatial.KDTree(x)
        self.max_d = 0

    def query(self, x, from2dto3d=False, N=32):
        if from2dto3d:
            d, indices = self.tree.query(x[0 : self.dim], N)
            self.max_d = max(np.max(d), self.max_d)
            U = np.zeros_like(self.U[0])
            for index in indices:
                U += self.U[index]
            return U / N
        else:
            try:
                d, indices = self.tree.query(x[0 : self.dim], N)
                self.max_d = max(d[0], self.max_d)
                rbf = RBFInterpolator(
                    self.x[indices],
                    self.U[indices],
                )
                return rbf(np.array([x[0 : self.dim]]))[0]
            except Exception as e:
                d, index = self.tree.query(x[0 : self.dim])
                self.max_d = max(d, self.max_d)
                return self.U[index]

    def max_distance(self):
        return self.max_d

    def reset(self):
        self.max_d = 0
