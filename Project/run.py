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


def read_data(t_idx_cut=0):
    # for name in names:
        # if name not in ['p', 'u', 'v', 'w']:
            # raise ValueError

    with open("./data/simu_params.txt", "r") as f:
        lines = f.readlines()
        params1 = [int(x) for x in lines[0].split(" ")]
        params2 = [float(x) for x in lines[1].split(" ")]
        params3 = [float(x) for x in lines[2].split(" ")]

    if t_idx_cut > 0:
        params1[0] = min(t_idx_cut, params1[0])
    
    nt, nx, ny = params1

    rows_1, rows_2 = (nt + 1) * nx * ny, (nt+1) * (nx+1) * (ny+1)

    p = np.loadtxt(f"./data/simu_p.txt", max_rows=rows_1).reshape((nt+1, nx, ny))
    u = np.loadtxt(f"./data/simu_u.txt", max_rows=rows_2).reshape((nt+1, nx+1, ny+1))
    v = np.loadtxt(f"./data/simu_v.txt", max_rows=rows_2).reshape((nt+1, nx+1, ny+1))
    w = np.loadtxt(f"./data/simu_w.txt", max_rows=rows_2).reshape((nt+1, nx+1, ny+1))
    
    # nt = 300//5
    # params1[0] = nt
    # f1 = f1[:nt+1]
    # f2 = f2[:nt+1]

    # if names[0] != 'p':
        # f1 = (f1[:, :-1, :-1] + f1[:, 1:, :-1] + f1[:, :-1, 1:] + f1[:, 1:, 1:]) / 4.
    # if names[1] != 'p':
        # f2 = (f2[:, :-1, :-1] + f2[:, 1:, :-1] + f2[:, :-1, 1:] + f2[:, 1:, 1:]) / 4.

    return params1, params2, params3, np.swapaxes(p, 1, 2), np.swapaxes(u, 1, 2), np.swapaxes(v, 1, 2), np.swapaxes(w, 1, 2)


def init():
    animated_items = []
    for ax, patch in zip(axs, patches):
        ax.add_patch(patch)
        animated_items.append(patch)

    pressure.set_data(p[0])
    animated_items.append(pressure)

    vorticity.set_data(w[0])
    animated_items.append(vorticity)

    if kwargs["quiver"]:
        uv.set_UVC(u_sampled[0], v_sampled[0], speed[0])
        animated_items.append(uv)

    time_text.set_text(time_template.format(0))
    animated_items.append(time_text)

    return animated_items


def update(t):
    animated_items = []

    for patch in patches[:-1]:
        patch.set_xy([din + delta_x[t], dbot + delta_y[t]])
        animated_items.append(patch)

    pressure.set_extent([delta_x[t], L+delta_x[t], 0 + delta_y[t], H + delta_y[t]])  # update position
    pressure.set_data(p[t])
    animated_items.append(pressure)

    vorticity.set_extent([delta_x[t], L+delta_x[t], 0 + delta_y[t], H + delta_y[t]])  # update position
    vorticity.set_data(w[t])
    animated_items.append(vorticity)

    if kwargs["quiver"]:
        uv.set_UVC(u_sampled[t], v_sampled[t], speed[t])
        animated_items.append(uv)

    time_text.set_text(time_template.format(t * save_dt))
    animated_items.append(time_text)
    return tuple(animated_items)

    # im.set_clim(vmin=np.amin(field[t]),vmax=np.amax(field[t]))


if __name__ == "__main__":

    kwargs = {"quiver":False}
    params1, params2, params3, p, u, v, w = read_data(t_idx_cut=0)

    nt, nx, ny = params1
    T, L, H, lbox, din, dbot, dt, h = params2
    t_start, alpha, strouhal, save_dt = params3
    n = int(nx // L)
    
    cmap1 = "Spectral_r"
    cmap2 = "Spectral_r"
    cmap3 = "viridis"
    cmap4 = "turbo_r"
    
    n_plots, win_h, trim, pad = (3, 9., 0.1, 0.05) if kwargs["quiver"] else (2, 8., 0., 0.05)
    
    t_array = np.linspace(0., T, nt + 1)

    # middle of cell : p value of MAC mesh
    xm = np.linspace(h/2., L-h/2., nx)
    ym = np.linspace(h/2., H-h/2., ny)
    
    # corners of cell : w value of MAC mesh
    x = np.linspace(0, L, nx + 1)
    y = np.linspace(0, H, ny + 1)

    xx_, yy_ = np.meshgrid(x, y)
    xxm, yym = np.meshgrid(xm, ym)

    # want u and v at corners of boxes
    # want w and p at center of boxes
    # w_full = np.sin(np.pi * xx_ / L) * np.sin(2 * np.pi * yy_ / H)
    # u = +1./H * np.sin(np.pi * xx_ / L) ** 2 * np.sin(2 * np.pi * yy_ / H)
    # v = -1./L * np.sin(2 * np.pi * xx_ / L) * np.sin(np.pi * yy_ / H) ** 2

    mask_middle = (din < xxm) * (xxm < din+lbox) * (dbot < yym) * (yym < dbot+1)
    mask_corner = (din < xx_) * (xx_ < din+lbox) * (dbot < yy_) * (yy_ < dbot+1)

    p_masked = np.ma.masked_array(p, [mask_middle for _ in range(nt+1)])
    w_masked = np.ma.masked_array(w, [mask_corner for _ in range(nt+1)])

    if kwargs["quiver"]:
        u = np.ma.masked_array(u, [mask_corner for _ in range(nt+1)])
        v = np.ma.masked_array(v, [mask_corner for _ in range(nt+1)])
        sample = max(1, n//5)
        x_sampled, y_sampled = xx_[::sample, ::sample], yy_[::sample, ::sample]
        u_sampled, v_sampled = u[:, ::sample, ::sample], v[:, ::sample, ::sample]
        speed = np.hypot(u_sampled, v_sampled)

    kappa = alpha / (2 * np.pi * strouhal)
    delta_x = 1. * kappa * (1. - np.cos(2. * np.pi * strouhal * t_array))
    delta_y = 0. * kappa * (1. - np.cos(2. * np.pi * strouhal * t_array))
    delta_x[t_array < t_start] = 0.
    delta_y[t_array < t_start] = 0.

    ################################  -  FIGURE  -  ################################
    fig, axs = plt.subplots(n_plots, 1, figsize=(12., 8.), constrained_layout=True, sharex="all", sharey="all")

    patches = [Rectangle((din, dbot), lbox, 1., fc="grey", ec="none", alpha=0.5, zorder=3) for _ in axs]

    pmin, pmax = np.amin(p), np.amax(p)
    pmin, pmax = -2., 2.
    wmin, wmax = np.amin(w)/10., np.amax(w)/10.
    # wmin, wmax = -15, 15

    pressure = axs[0].imshow(p[0], extent=(0, L, 0, H), vmin=pmin, vmax=pmax, cmap=cmap1, origin="lower")
    vorticity = axs[1].imshow(w[0], extent=(0, L, 0, H), vmin=wmin, vmax=wmax, cmap=cmap2, origin="lower")
    cbar_p = fig.colorbar(pressure, ax=axs[0], pad=pad)
    cbar_w = fig.colorbar(vorticity, ax=axs[1], pad=pad)

    if kwargs["quiver"]:
        uv = axs[2].quiver(x_sampled, y_sampled, u_sampled[0], v_sampled[0], speed[0], cmap=cmap4, width=0.0015, zorder=2) # vmin=np.amin(speed), vmax=np.amax(speed)  
        cbar_q = fig.colorbar(uv, ax=axs[2], pad=pad)
        # stream = ax.streamplot(x, y, u, v, density=1.5)
    
    for ax in axs:
        ax.axis([0., L, 0. + trim*H, H - trim*H])
        ax.set_aspect('equal')
    
    time_template = "t = {:5.3f} [s]"
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
