import numpy as np


def walldistance(coord):
    # the ramp starts at x = 0
    x, y, _ = coord
    # y = 0
    d1 = y
    # y = tan(24) * x
    # tan(24) * x - y = 0
    k2 = np.tan(24 * np.pi / 180)
    d2 = np.abs(k2 * x - y) / np.hypot(k2, 1)
    return min(d1, d2)
