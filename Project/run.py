import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse as sp
from matplotlib.patches import Rectangle
from matplotlib.animation import FuncAnimation

ftSz1, ftSz2, ftSz3 = 20, 17, 15
plt.rcParams['font.family'] = 'monospace'


def read_data_neumann_example(filename):
    values = np.loadtxt(filename)
    return values.T


def read_data():
    with open("./data/simu_params.txt", "r") as f:
        lines = f.readlines()
        params1 = [int(x) for x in lines[0].split(" ")]
        params2 = [float(x) for x in lines[1].split(" ")]
        # param_list = (nt, nx, ny, T, L, H, lbox, din, dbot, dt, h, alpha, strouhal)
    
    nt, nx, ny = params1
    u = np.loadtxt("./data/simu_u.txt").reshape((nt+1, nx+1, ny+1))
    v = np.loadtxt("./data/simu_v.txt").reshape((nt+1, nx+1, ny+1))
    w = np.loadtxt("./data/simu_w.txt").reshape((nt+1, nx+1, ny+1))
    p = np.loadtxt("./data/simu_p.txt").reshape((nt+1, nx+0, ny+0))

    return params1, params2, np.swapaxes(u, 1, 2), np.swapaxes(v, 1, 2), np.swapaxes(w, 1, 2), np.swapaxes(p, 1, 2)


def init():
    pressure.set_data(p[0])
    vorticity.set_data(w[0])
    uv.set_UVC(u[0], v[0])

    time_text.set_text(time_template.format(0))
    return pressure, vorticity, uv, time_text


def update(t):
    pressure.set_data(p[t])
    vorticity.set_data(w[t])
    uv.set_UVC(u[t], v[t])

    time_text.set_text(time_template.format(t * dt))
    return pressure, vorticity, uv, time_text

    # im.set_clim(vmin=np.amin(field[t]),vmax=np.amax(field[t]))


if __name__ == "__main__":

    params1, params2, u, v, w, p = read_data()
    nt, nx, ny = params1
    T, L, H, lbox, din, dbot, dt, h, alpha, strouhal = params2
    n = int(nx // L)
    
    # L, H = 15, 5
    # n = 5
    # h = 1. / n

    # middle of cell : p value of MAC mesh
    xm = np.linspace(h/2., L-h/2., nx)
    ym = np.linspace(h/2., H-h/2., ny)
    
    # corners of cell : w value of MAC mesh
    x = np.linspace(0, L, nx + 1)
    y = np.linspace(0, H, ny + 1)

    xx_, yy_ = np.meshgrid(x, y)
    xx, yy = np.meshgrid(xm, ym)

    # want u and v at corners of boxes
    # want w and p at center of boxes

    # w_full = np.sin(np.pi * xx_ / L) * np.sin(2 * np.pi * yy_ / H)
    # u = +1./H * np.sin(np.pi * xx_ / L) ** 2 * np.sin(2 * np.pi * yy_ / H)
    # v = -1./L * np.sin(2 * np.pi * xx_ / L) * np.sin(np.pi * yy_ / H) ** 2

    w = (w[:, :-1, :-1] + w[:, 1:, :-1] + w[:, :-1, 1:] + w[:, 1:, 1:]) / 4.

    mask = (din < xx) * (xx < din+lbox) * (dbot < yy) * (yy < dbot+1)
    p = np.ma.masked_array(p, [mask for _ in range(nt+1)])
    w = np.ma.masked_array(w, [mask for _ in range(nt+1)])

    mask_ = (din < xx_) * (xx_ < din+lbox) * (dbot < yy_) * (yy_ < dbot+1)
    u = np.ma.masked_array(u, [mask_ for _ in range(nt+1)])
    v = np.ma.masked_array(v, [mask_ for _ in range(nt+1)])
    # field[:, mask] = 0.


    ################################  -  FIGURE  -  ################################
    fig, axs = plt.subplots(2, 1, figsize=(12., 8.), constrained_layout=True)

    for ax in axs:
        rect = Rectangle((3., 2.), 5., 1., fc="grey", ec="none")
        ax.add_patch(rect)

    # ax.pcolormesh(x, y, u, cmap="bwr", shading="flat", aspect="")

    pmin, pmax = np.amin(p), np.amax(p)    
    pressure = axs[0].imshow(w[0], extent=(0, L, 0, H), vmin=pmin, vmax=pmax, cmap="Spectral", origin="lower")
    cbar_p = plt.colorbar(pressure, ax=axs[0])

    wmin, wmax = np.amin(w), np.amax(w)
    vorticity = axs[1].imshow(w[0], extent=(0, L, 0, H), vmin=wmin, vmax=wmax, cmap="Spectral", origin="lower")
    cbar_w = plt.colorbar(vorticity, ax=axs[1])

    uv = ax.quiver(xx_, yy_, u[0], v[0], width=0.0015, zorder=2)


    # stream = ax.streamplot(x, y, u, v, density=1.5)

    # im.set_data(field)  # update values
    # im.set_extent([-1, L-1, 0, H])  # update position

    # ax.axis([0., L, 0., H])
    
    time_template = "t = {:4.2f} [s]"
    time_text = axs[0].text(0.8,0.9, time_template.format(0),
                        fontsize=ftSz2, transform=axs[0].transAxes,
                        bbox=dict(boxstyle="round", fc="wheat", ec="none", alpha=0.85))

    anim = FuncAnimation(fig, update, nt+1, interval=500, blit=False, init_func=init, repeat_delay=3000)

    plt.show()
