def walldistance(coord, _):
    # step height is 0.0375, step corner is at x = 0, y = 0
    h = 0.0375
    x, y, _ = coord
    if x <= 0:
        return min(4 * h - y, y)
    elif y >= 0:
        d1 = 4 * h - y
        d2 = y + h
        d3 = (x**2 + y**2) ** 0.5
        return min(d1, d2, d3)
    else:
        d1 = y + h
        d2 = x
        return min(d1, d2)


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    import numpy as np

    h = 0.0375
    x = np.linspace(-4 * h, 20 * h, 201)
    y = np.linspace(-h, 4 * h, 101)
    x, y = np.meshgrid(x, y)
    d = np.vectorize(lambda x, y: walldistance([x, y, 0]))(x, y)
    for i in range(len(x)):
        for j in range(len(y)):
            if x[i, j] < 0 and y[i, j] < 0:
                d[i, j] = np.nan

    fig, ax = plt.subplots(figsize=(6, 3), layout="constrained")
    c = ax.contourf(x / h, y / h, d / h)
    fig.colorbar(c)

    ax.set_aspect("equal")

    ax.tick_params(which="both", direction="in")
    ax.set_xlabel(r"$x / h$")
    ax.set_ylabel(r"$y / h$")

    fig.savefig("bfs_walldistance.png", dpi=300)
