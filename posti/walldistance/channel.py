def walldistance(coord):
    # half height is 1, centerline at y = 0
    _, y, _ = coord
    if y >= 0:
        return 1 - y
    else:
        return y + 1