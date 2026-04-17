import numpy as np


def walldistance(coord, args):
    if "theta" in args:
        theta = args["theta"]
    else:
        theta = 24
    x, y, _ = coord
    d1 = y
    k2 = np.tan(theta * np.pi / 180)
    d2 = np.abs(k2 * x - y) / np.hypot(k2, 1)
    return min(d1, d2)
