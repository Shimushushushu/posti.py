def walldistance(coord, _):
    x, y, _ = coord
    if x >= 0:
        return y
    else:
        return (x**2 + y**2)**0.5