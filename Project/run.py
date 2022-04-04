import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse as sp
from time import perf_counter
from matplotlib.patches import Rectangle
from matplotlib.animation import FuncAnimation

ftSz1, ftSz2, ftSz3 = 20, 17, 15
plt.rcParams['font.family'] = 'monospace'

path_dir = "./data"
size_p, size_o = None, None

def read_data_neumann_example(filename):
    values = np.loadtxt(filename)
    return values.T


def read_data_init():

    with open("./data/simu_params.txt", "r") as f:
        lines = f.readlines()
        params1 = [int(x) for x in lines[0].split(" ")]
        params2 = [float(x) for x in lines[1].split(" ")]
        params3 = [float(x) for x in lines[2].split(" ")]
    
    nt, nx, ny, n, save_modulo = params1
    params2[1] *= save_modulo  # dt multiplied since not all times are saved

    global size_p, size_o
    size_p, size_o = nx * ny, (nx+1) * (ny+1)

    # p = np.loadtxt(f"./data/simu_p.txt", max_rows=(nt+1) * size_p).reshape((nt+1, nx, ny))
    # u = np.loadtxt(f"./data/simu_u.txt", max_rows=(nt+1) * size_o).reshape((nt+1, nx+1, ny+1))
    # v = np.loadtxt(f"./data/simu_v.txt", max_rows=(nt+1) * size_o).reshape((nt+1, nx+1, ny+1))
    # w = np.loadtxt(f"./data/simu_w.txt", max_rows=(nt+1) * size_o).reshape((nt+1, nx+1, ny+1))
    
    return params1, params2, params3


def read_block(idx_start, n_steps):
    exec_start = perf_counter()
    print("Loading from  {:4d}  to  {:5d}".format(idx_start, idx_start+n_steps), end="\r")
    n_steps = min(n_steps, (nt + 1) - idx_start)
    skip, maxr = idx_start * size_p, n_steps * size_p
    p = np.loadtxt(f"{path_dir}/simu_p.txt", skiprows=skip, max_rows=maxr).reshape((n_steps, nx, ny))
    
    skip, maxr = idx_start * size_o, n_steps * size_o
    u = np.loadtxt(f"{path_dir}/simu_u.txt", skiprows=skip, max_rows=maxr).reshape((n_steps, nx+1, ny+1))
    v = np.loadtxt(f"{path_dir}/simu_v.txt", skiprows=skip, max_rows=maxr).reshape((n_steps, nx+1, ny+1))
    w = np.loadtxt(f"{path_dir}/simu_w.txt", skiprows=skip, max_rows=maxr).reshape((n_steps, nx+1, ny+1))

    p, w = np.swapaxes(p, 1, 2), np.swapaxes(w, 1, 2)
    u, v = np.swapaxes(u, 1, 2), np.swapaxes(v, 1, 2)

    p_masked = np.ma.masked_array(p, [mask_middle for _ in p])
    w_masked = np.ma.masked_array(w, [mask_corner for _ in w])

    # if kwargs["quiver"]:
    #     u = np.ma.masked_array(u, [mask_corner for _ in u])
    #     v = np.ma.masked_array(v, [mask_corner for _ in v])
    #     sample = max(1, n//5)
    #     x_sampled, y_sampled = xx_[::sample, ::sample], yy_[::sample, ::sample]
    #     u_sampled, v_sampled = u[:, ::sample, ::sample], v[:, ::sample, ::sample]
    #     speed = np.hypot(u_sampled, v_sampled)
    print(f"Loading from  {idx_start:4d}  to  {idx_start+n_steps:5d}  in {perf_counter() - exec_start:.3f} s")

    return p, u, v, w, p_masked, w_masked


def find_bounds(filename):
    min_value, max_value = np.inf, -np.inf
    with open(filename, "r") as f:
        for line in f:
            value = float(line)
            min_value = min(min_value, value)
            max_value = max(max_value, value)
    return min_value, max_value


def find_all_bounds():
    bounds_p = find_bounds(f"{path_dir}/simu_p.txt")
    bounds_u = find_bounds(f"{path_dir}/simu_u.txt")
    bounds_v = find_bounds(f"{path_dir}/simu_v.txt")
    bounds_w = find_bounds(f"{path_dir}/simu_w.txt")
    return bounds_p, bounds_u, bounds_v, bounds_w


def init():
    animated_items = []
    for ax, patch in zip(axs, patches):
        ax.add_patch(patch)
        animated_items.append(patch)

    pressure.set_data(p[0])
    animated_items.append(pressure)

    vorticity.set_data(w[0])
    animated_items.append(vorticity)

    # if kwargs["quiver"]:
    #     uv.set_UVC(u_sampled[0], v_sampled[0], speed[0])
    #     animated_items.append(uv)

    time_text.set_text(time_template.format(0))
    animated_items.append(time_text)

    return animated_items


def update(t):
    global current_t_idx, p, u, v, w, p_masked, w_masked

    if (current_t_idx > 0) and (t == 0):
        current_t_idx = 0
        p, u, v, w, p_masked, w_masked = read_block(current_t_idx, block_size)

    if current_t_idx + block_size <= t:
        current_t_idx += block_size
        p, u, v, w, p_masked, w_masked = read_block(current_t_idx, block_size)

    animated_items = []

    for patch in patches:
        patch.set_xy([din + delta_x[t], dbot + delta_y[t]])
        animated_items.append(patch)
     
    pressure.set_extent([delta_x[t], L+delta_x[t], 0 + delta_y[t], H + delta_y[t]])  # update position
    pressure.set_data(p_masked[t % block_size])
    animated_items.append(pressure)

    vorticity.set_extent([delta_x[t], L+delta_x[t], 0 + delta_y[t], H + delta_y[t]])  # update position
    vorticity.set_data(w_masked[t % block_size])
    animated_items.append(vorticity)

    # if kwargs["quiver"]:
    #     uv.set_UVC(u_sampled[t], v_sampled[t], speed[t])
    #     animated_items.append(uv)

    time_text.set_text(time_template.format(t * dt))
    animated_items.append(time_text)
    return tuple(animated_items)

    # im.set_clim(vmin=np.amin(field[t]),vmax=np.amax(field[t]))


if __name__ == "__main__":

    kwargs = {"quiver":False}
    block_size = 20

    ############################################################################################
    ####################################  -  SETUP DATA  -  ####################################
    
    params1, params2, params3 = read_data_init()
    # t1 = perf_counter()
    # (pmin, pmax), (umin, umax), (vmin, vmax), (wmin, wmax) = find_all_bounds()
    # print(f"Time to find bound = {perf_counter() - t1:.3f} s")

    nt, nx, ny, n, _ = params1
    T, dt, h, L, H, lbox, din, dbot = params2
    swing_start, pert_start, pert_dt, alpha, strouhal = params3
    
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

    xx_, yy_ = np.meshgrid(x, y)        # want u and v at corners of boxes
    xxm, yym = np.meshgrid(xm, ym)      # want w and p at center of boxes

    mask_middle = (din < xxm) * (xxm < din+lbox) * (dbot < yym) * (yym < dbot+1)
    mask_corner = (din <= xx_) * (xx_ <= din+lbox) * (dbot <= yy_) * (yy_ <= dbot+1.)

    kappa = alpha / (2 * np.pi * strouhal)
    delta_x = 0. * kappa * (1. - np.cos(2. * np.pi * strouhal * t_array))
    delta_y = 0.5 * pert_dt / (2. * np.pi) * (1. - np.cos(2. * np.pi * (t_array - pert_start) / pert_dt))

    delta_x[t_array < swing_start] = 0.
    delta_y[~((pert_start < t_array) * (t_array < pert_start + pert_dt))] = 0.

    current_t_idx = 0
    p, u, v, w, p_masked, w_masked = read_block(current_t_idx, block_size)


    #############################################################################################
    ######################################  -  FIGURE  -  #######################################
    fig, axs = plt.subplots(n_plots, 1, figsize=(12., 8.), constrained_layout=True, sharex="all", sharey="all")

    # pmin, pmax = np.amin(p), np.amax(p)
    pmin, pmax = -2., 2.
    # wmin, wmax = np.amin(w)/10., np.amax(w)/10.
    wmin, wmax = -25, 25
    # pmin, pmax = np.array([pmin, pmax]) / 1.
    # wmin, wmax = np.array([wmin, wmax]) / 5.

    pressure = axs[0].imshow(p[0], extent=(0, L, 0, H), vmin=pmin, vmax=pmax, cmap=cmap1, origin="lower")
    vorticity = axs[1].imshow(w[0], extent=(0, L, 0, H), vmin=wmin, vmax=wmax, cmap=cmap2, origin="lower")
    cbar_p = fig.colorbar(pressure, ax=axs[0], pad=pad)
    cbar_w = fig.colorbar(vorticity, ax=axs[1], pad=pad)

    patches = [Rectangle((din, dbot), lbox, 1., fc="none", ec="black", alpha=1., zorder=1) for _ in axs]

    # if kwargs["quiver"]:
    #     uv = axs[2].quiver(x_sampled, y_sampled, u_sampled[0], v_sampled[0], speed[0], cmap=cmap4, width=0.0015, zorder=2) # vmin=np.amin(speed), vmax=np.amax(speed)  
    #     cbar_q = fig.colorbar(uv, ax=axs[2], pad=pad)
    #     # stream = ax.streamplot(x, y, u, v, density=1.5)
    
    for ax in axs:
        ax.axis([0., L, 0. + trim*H, H - trim*H])
        # ax.axis([2., 10., 1.5, 3.5])
        ax.set_aspect('equal')
    
    time_template = "t = {:5.3f} [s]"
    time_text = axs[0].text(0.8,0.9, time_template.format(0),
                        fontsize=ftSz2, transform=axs[0].transAxes,
                        bbox=dict(boxstyle="round", fc="wheat", ec="none", alpha=0.85))

    anim_type = 0

    if anim_type == 0:
        anim = FuncAnimation(fig, update, nt+1, interval=1000, blit=True, init_func=init, repeat=True, repeat_delay=2000)
    elif anim_type == 1:
        update(nt)
        # stream = axs[0].streamplot(x, y, u[-1], v[-1], density=1)
    elif anim_type == 3:
        init()
        for t in range(nt+1):
            update(t)
            fig.savefig("./anim/frame_{:05d}.png", format="png", bbox_inches='tight')


    plt.show()
