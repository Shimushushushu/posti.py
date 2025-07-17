import numpy as np
from scipy.optimize import minimize_scalar


def faith_geometry(r):
    if r <= 3 / 2:
        return (np.cos((2 / 3) * np.pi * r) + 1) / 2
    else:
        return 0


def faith_distance_squared(t, r, y):
    return (faith_geometry(t) - y) ** 2 + (t - r) ** 2


def walldistance(coord):
    # hill is located at x = 0, z = 0, hill height = 1
    x, y, z = coord
    r = (x**2 + z**2) ** 0.5
    t = minimize_scalar(
        lambda t: faith_distance_squared(t, r, y), bounds=(0, r), method="bounded"
    ).x
    d1 = faith_distance_squared(t, r, y)**0.5
    d2 = 5.3 - y
    d3 = min(4 - z, z + 4)
    return min(d1, d2, d3)
