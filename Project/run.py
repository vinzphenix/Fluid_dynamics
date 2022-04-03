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


def read_data(names):
    for name in names:
        if name not in ['p', 'u', 'v', 'w']:
            raise ValueError

    with open("./data/simu_params.txt", "r") as f:
        lines = f.readlines()
        params1 = [int(x) for x in lines[0].split(" ")]
        params2 = [float(x) for x in lines[1].split(" ")]
        # param_list = (nt, nx, ny, T, L, H, lbox, din, dbot, dt, h, alpha, strouhal)
    
    nt, nx, ny = params1
    inc = [int(name != 'p') for name in names]

    f1 = np.loadtxt(f"./data/simu_{names[0]}.txt").reshape((nt+1, nx+inc[0], ny+inc[0]))
    f2 = np.loadtxt(f"./data/simu_{names[1]}.txt").reshape((nt+1, nx+inc[1], ny+inc[1]))
    
    # nt = 300//5
    # params1[0] = nt
    # f1 = f1[:nt+1]
    # f2 = f2[:nt+1]

    # if names[0] != 'p':
        # f1 = (f1[:, :-1, :-1] + f1[:, 1:, :-1] + f1[:, :-1, 1:] + f1[:, 1:, 1:]) / 4.
    # if names[1] != 'p':
        # f2 = (f2[:, :-1, :-1] + f2[:, 1:, :-1] + f2[:, :-1, 1:] + f2[:, 1:, 1:]) / 4.

    return params1, params2, np.swapaxes(f1, 1, 2), np.swapaxes(f2, 1, 2)


def init():
    pressure.set_data(f1[0])
    vorticity.set_data(f2[0])
    # uv.set_UVC(u[0], v[0])

    time_text.set_text(time_template.format(0))
    return pressure, vorticity, time_text, # uv


def update(t):
    pressure.set_data(f1[t])
    vorticity.set_data(f2[t])
    # uv.set_UVC(u[t], v[t])

    time_text.set_text(time_template.format(t * save_dt))
    return pressure, vorticity, time_text, # uv

    # im.set_clim(vmin=np.amin(field[t]),vmax=np.amax(field[t]))


if __name__ == "__main__":

    params1, params2, f1, f2 = read_data(['p', 'w'])

    nt, nx, ny = params1
    T, L, H, lbox, din, dbot, dt, h, alpha, strouhal, save_dt = params2
    n = int(nx // L)
    
    cmap1 = "Spectral_r"
    cmap2 = "Spectral_r"
    
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

    mask = (din < xx) * (xx < din+lbox) * (dbot < yy) * (yy < dbot+1)
    mask_ = (din < xx_) * (xx_ < din+lbox) * (dbot < yy_) * (yy_ < dbot+1)

    f1_masked = np.ma.masked_array(f1, [mask for _ in range(nt+1)])
    f2_masked = np.ma.masked_array(f2, [mask_ for _ in range(nt+1)])


    ################################  -  FIGURE  -  ################################
    fig, axs = plt.subplots(2, 1, figsize=(12., 8.), constrained_layout=True, sharex="all", sharey="all")

    for ax in axs:
        rect = Rectangle((3., 2.), 5., 1., fc="grey", ec="none", alpha=0.5, zorder=3)
        ax.add_patch(rect)


    f1min, f1max = np.amin(f1), np.amax(f1)
    f1min, f1max = -2., 2.
    pressure = axs[0].imshow(f1[0], extent=(0, L, 0, H), vmin=f1min, vmax=f1max, cmap=cmap1, origin="lower")
    cbar_p = plt.colorbar(pressure, ax=axs[0])

    # f2min, f2max = np.amin(f2)/1., np.amax(f2)/1.
    f2min, f2max = -15, 15
    vorticity = axs[1].imshow(f2[0], extent=(0, L, 0, H), vmin=f2min, vmax=f2max, cmap=cmap2, origin="lower")
    cbar_w = plt.colorbar(vorticity, ax=axs[1])

    # uv = ax.quiver(xx_, yy_, u[0], v[0], width=0.0015, zorder=2)
    # stream = ax.streamplot(x, y, u, v, density=1.5)

    # im.set_data(field)  # update values
    # im.set_extent([-1, L-1, 0, H])  # update position
    # ax.axis([0., L, 0., H])
    
    time_template = "t = {:6.4f} [s]"
    time_text = axs[0].text(0.8,0.9, time_template.format(0),
                        fontsize=ftSz2, transform=axs[0].transAxes,
                        bbox=dict(boxstyle="round", fc="wheat", ec="none", alpha=0.85))

    anim_type = 0

    if anim_type == 0:
        anim = FuncAnimation(fig, update, nt+1, interval=500, blit=False, init_func=init, repeat=False)
    elif anim_type == 1:
        update(nt)
        # stream = axs[0].streamplot(x, y, u[-1], v[-1], density=1)
    elif anim_type == 3:
        init()
        for t in range(nt+1):
            update(t)
            fig.savefig("./anim/frame_{:05d}.png", format="png", bbox_inches='tight')


    plt.show()
