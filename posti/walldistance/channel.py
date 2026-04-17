def walldistance(coord, args):
    if "h" in args:
        h = args["h"]
    else:
        h = 1
    _, y, _ = coord
    return min(h - y, y + h)
