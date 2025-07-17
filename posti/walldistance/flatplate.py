def walldistance(coord):
    # flatplate starts at x = 0
    x, y, _ = coord
    if x >= 0:
        return y
    else:
        return (x**2 + y**2)**0.5