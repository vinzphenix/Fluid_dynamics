import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

ftSz1, ftSz2, ftSz3 = 20, 17, 15
plt.rcParams['font.family'] = 'monospace'


if __name__ == "__main__":
    L, H = 15, 5
    n = 10
    h = 1. / n

    # middle of cell : p value of MAC mesh
    xm = np.linspace(h/2., L-h/2., L * n)
    ym = np.linspace(h/2., H-h/2., H * n)
    
    # corners of cell : w value of MAC mesh
    x = np.linspace(0, L, L*n + 1)
    y = np.linspace(0, H, H*n + 1)

    xx_, yy_ = np.meshgrid(x, y)
    xx, yy = np.meshgrid(xm, ym)

    # want u and v at corners of boxes
    # want w and p at center of boxes
    w = np.sin(np.pi * xx / L) * np.sin(2 * np.pi * yy / H)
    u = +1./H * np.sin(np.pi * xx_ / L) ** 2 * np.sin(2 * np.pi * yy_ / H)
    v = -1./L * np.sin(2 * np.pi * xx_ / L) * np.sin(np.pi * yy_ / H) ** 2
    
    mask = (3 < xx) * (xx < 8) * (2 < yy) * (yy < 3)
    w[mask] = 0.5

    fig, ax = plt.subplots(1, 1, figsize=(12., 6.), constrained_layout=True)

    rect = Rectangle((3., 2.), 5., 1., fc="grey", ec="none", zorder=3)
    ax.add_patch(rect)

    # ax.pcolormesh(x, y, u, cmap="bwr", shading="flat", aspect="")
    vorticity = ax.imshow(w, extent=(0, L, 0, H), cmap="bwr", origin="lower")
    stream = ax.streamplot(x, y, u, v, density=1.5)

    # vorticity.set_data(w)  # update values
    # vorticity.set_extent([-1, L-1, 0, H])  # update position

    # ax.axis([0., L, 0., H])

    plt.show()
