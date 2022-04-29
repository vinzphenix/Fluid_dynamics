import matplotlib.pyplot as plt
import numpy as np

ftSz1, ftSz2, ftSz3 = 20, 17, 15
plt.rcParams['font.family'] = 'monospace'
plt.rcParams["text.usetex"] = False

if __name__ == "__main__":
    fig, ax = plt.subplots(1, 1, figsize=(10., 6.))
    fig.tight_layout()

    h = 0.1
    nx, ny = 7, 3
    L, H = nx * h, ny * h

    x_w, y_w = np.linspace(0., L, nx + 1), np.linspace(0., H, ny + 1)
    x_u, y_u = x_w, y_w - h / 2.
    x_v, y_v = x_w - h / 2., y_w
    x_t, y_t = x_w - h / 2., y_w - h / 2.
    x_grid, y_grid = np.r_[-h / 2., x_w, x_w[-1] + h / 2.], np.r_[-h / 2., y_w, y_w[-1] + h / 2.]

    x_w, y_w = np.meshgrid(x_w, y_w)
    x_u, y_u = np.meshgrid(x_u, y_u)
    x_v, y_v = np.meshgrid(x_v, y_v)
    x_t, y_t = np.meshgrid(x_t, y_t)
    x_grid, y_grid = np.meshgrid(x_grid, y_grid)

    ax.plot(x_grid[:, 1:-1], y_grid[:, 1:-1], ls=':', color='grey')
    ax.plot(x_grid.T[:, 1:-1], y_grid.T[:, 1:-1], ls=':', color='grey')
    ax.plot(x_w[:, 0], y_w[:, 0], ls='-', color="black", lw=3)  # left wall
    ax.plot(x_w[0, :], y_w[0, :], ls='-', color="black", lw=3)  # bottom wall

    ms = 10
    ax.plot(x_u[1:, :], y_u[1:, :], ">", markersize=ms, color='C0', label=r'$u$')
    ax.plot(x_v[:, 1:], y_v[:, 1:], "^", markersize=ms, color='C2', label=r'$v$')
    ax.plot(x_w, y_w, "o", markersize=ms, color='black', label=r'$\omega$')
    ax.plot(x_t[1:, 1:], y_t[1:, 1:], "*", markersize=ms, color='C1', label=r'$T$')
    ax.scatter(x_t[1:, 1:], y_t[1:, 1:], s=150, fc='none', ec='red', label=r'$p$')

    # ghosts
    alpha = 0.35
    ax.plot(x_u[0, :], y_u[0, :], ">", markersize=ms, color='C0', alpha=alpha)
    ax.plot(x_v[:, 0], y_v[:, 0], "^", markersize=ms, color='C2', alpha=alpha)
    ax.plot(x_t[0, :], y_t[0, :], "*", markersize=ms, color='C1', alpha=alpha)
    ax.plot(x_t[1:, 0], y_t[1:, 0], "*", markersize=ms, color='C1', alpha=alpha)

    dots_horizontal = np.linspace(L + h, L + 1.5 * h, 3)
    dots_vertical = np.linspace(H + 0.75 * h, H + 1.25 * h, 3)
    ax.plot(dots_horizontal, H * 0.6 * np.ones_like(dots_horizontal), "o", markersize=4, color="grey")
    ax.plot(L * 0.5 * np.ones_like(dots_vertical), dots_vertical, "o", markersize=4, color="grey")
    ax.plot(dots_horizontal, dots_vertical, "o", markersize=4, color="grey")

    lines_labels = [ax.get_legend_handles_labels()]
    lines, labels = [sum(lol, []) for lol in zip(*lines_labels)]
    lines, labels = lines[2::nx], labels[2::nx]

    ax.legend(lines, labels, loc='lower right', fontsize=ftSz1)
    ax.set_aspect("equal")
    ax.axis("off")

    fig_name = "mac_mesh"
    fig.savefig(f"../figures/{fig_name}.svg", format="svg", bbox_inches='tight')

    plt.show()
