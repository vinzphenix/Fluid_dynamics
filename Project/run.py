import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from tqdm import tqdm
from scipy.integrate import cumulative_trapezoid
from scipy.interpolate import RegularGridInterpolator
from time import perf_counter
from matplotlib.animation import FuncAnimation, PillowWriter, FFMpegWriter
from matplotlib.patches import Rectangle
from mpl_toolkits.axes_grid1 import make_axes_locatable


ftSz1, ftSz2, ftSz3 = 20, 17, 15
plt.rcParams['font.family'] = 'monospace'

path_dir = "./data"
global size_p, size_u, size_v


def read_data_init():

    with open("./data/simu_params.txt", "r") as f:
        lines = f.readlines()
        params1 = [int(x) for x in lines[0].split(" ")]
        params2 = [float(x) for x in lines[1].split(" ")]
        params3 = [float(x) for x in lines[2].split(" ")]
    
    nt, nx, ny, n, save_modulo = params1
    params2[1] *= save_modulo  # dt multiplied since not all times are saved

    global size_p, size_u, size_v
    size_p, size_u, size_v = nx * ny, (nx+1) * (ny+2), (nx+2) * (ny+1)

    file_ptrs = [open(f"{path_dir}/simu_p.txt"), open(f"{path_dir}/simu_u.txt"), open(f"{path_dir}/simu_v.txt")]
    
    df_p = pd.read_csv("./data/simu_p.txt", chunksize=size_p, header=None)
    df_u = pd.read_csv("./data/simu_u.txt", chunksize=size_u, header=None)
    df_v = pd.read_csv("./data/simu_v.txt", chunksize=size_v, header=None)
    
    # find bounds
    # t1 = perf_counter()
    # (pmin, pmax), (umin, umax), (vmin, vmax), (wmin, wmax) = find_all_bounds()
    # print(f"Time to find bound = {perf_counter() - t1:.3f} s")

    return params1, params2, params3, [df_p, df_u, df_v]


def read_block(dataframe_list):
    df_p, df_u, df_v = dataframe_list
    
    p = np.array(next(df_p)).reshape((nx, ny))
    u = np.array(next(df_u)).reshape((nx+1, ny+2))
    v = np.array(next(df_v)).reshape((nx+2, ny+1))
    
    w = ((v[1:, :] - v[:-1, :]) - (u[:, 1:] - u[:, :-1])) / h  # dv/dx - du/dy --> values at corners
    w = 0.25 * (w[1:, 1:] + w[1:, :-1] + w[:-1, 1:] + w[:-1, :-1])  # average w on corners to get value in the middle

    p_masked = np.ma.masked_array(p, mask_middle)
    w_masked = np.ma.masked_array(w, mask_middle)

    psi = np.zeros((nx+1, ny+1))
    psi[:, :ny//2+1] = cumulative_trapezoid(u[:, 1:ny//2+2], dx=h, initial=0., axis=1)
    psi[:, ny:ny//2:-1] = cumulative_trapezoid(u[:, ny-1:ny//2-1:-1], dx=h, initial=0., axis=1)
    psi = np.ma.masked_array(psi.T, mask_corner)
    # psi[:, :ny//2] = h * np.cumsum(u[:, 1:ny//2+1], axis=1)  # more readable and "rectangle" integration
    # for j in range(ny//2):                                   # equivalent with for loops
        # psi[:, j+1] = psi[:, j] + h * u[:, j+1]
    # for j in range(ny, ny//2, -1):
        # psi[:, j-1] = psi[:, j] + h * u[:, j-1]

    return p.T, u, v, w.T, p_masked.T, w_masked.T, psi


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

    animated_items.append(pressure)
    animated_items.append(vorticity)

    if kwargs["stream"]:
        animated_items.append(strm[0].collections)
    if kwargs["streak"]:
        for j in range(nStreakLines):
            streaklines_x[j] = -np.ones(nParticles)
            streaklines_y[j] = spots_y[j]*np.ones(nParticles)        
        for j, line in enumerate(lines):
            line.set_data([], [])
        animated_items += lines

    time_text.set_text(time_str.format(0))
    animated_items.append(time_text)

    return animated_items


def update(t):
    global p, u, v, w, p_masked, w_masked, psi

    if (t > 0):
        p, u, v, w, p_masked, w_masked, psi = read_block(dataframes)

    animated_items = []

    for patch in patches:
        patch.set_xy([din + delta_x[t], dbot + delta_y[t]])
        animated_items.append(patch)
     
    pressure.set_extent([delta_x[t], L+delta_x[t], 0 + delta_y[t], H + delta_y[t]])  # update position
    pressure.set_data(p)
    animated_items.append(pressure)

    vorticity.set_extent([delta_x[t], L+delta_x[t], 0 + delta_y[t], H + delta_y[t]])  # update position
    vorticity.set_data(w)
    animated_items.append(vorticity)

    if kwargs["stream"]:
        for item in strm[0].collections:
            item.remove()
        strm[0] = axs[0].contour(xx_, yy_, psi, colors=strm_color, alpha=strm_alpha, levels=strm_levels)
        animated_items.append(strm[0].collections)
    
    if kwargs["streak"]:
        U, V = 0.5 * (u[:, :-1] + u[:, 1:]), 0.5 * (v[:-1, :] + v[1:, :])
        u_interp = RegularGridInterpolator((x, y), U, bounds_error=False, fill_value=0)
        v_interp = RegularGridInterpolator((x, y), V, bounds_error=False, fill_value=0)
        cycle = int(np.ceil((L / 1.) / nParticles / dt)) # [temps pour traverser / nbre particules] = delai
        
        for j, (line, x_particles, y_particles) in enumerate(zip(lines, streaklines_x, streaklines_y)):
            if t % cycle == 0:
                x_particles[(t//cycle) % nParticles] = 0.
                y_particles[(t//cycle) % nParticles] = spots_y[j]
            mask = (0. <= x_particles) * (x_particles <= L) * (0. <= y_particles) * (y_particles <= H)
            args = np.c_[x_particles[mask], y_particles[mask]]
            
            x_particles[mask] += dt * u_interp(args)  # Euler explicite
            y_particles[mask] += dt * v_interp(args)  # Euler explicite

            x_particles[~mask] = np.nan
            y_particles[~mask] = np.nan
            line.set_data(x_particles, y_particles)
        
        animated_items += lines


    time_text.set_text(time_str.format(t * dt))
    animated_items.append(time_text)
    return tuple(animated_items)

    # im.set_clim(vmin=np.amin(field[t]),vmax=np.amax(field[t]))


def make_colorbar_with_padding(ax):
    """
    Create colorbar axis that fits the size of a plot - detailed here: http://chris35wills.github.io/matplotlib_axis/
    """
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="2%", pad=0.1)
    return(cax) 


def modify_html(html_filename):
    with open(html_filename, "r") as f:
        html_lines = f.readlines()
        old_line = html_lines[185]
        new_line = old_line[:15] + str(nt+1) + ";\n"
        html_lines[185] = new_line

    with open(html_filename, "w") as f:
        f.write("".join(html_lines))
        f.close()




if __name__ == "__main__":

    kwargs = {"stream":True, "streak":True}

    ############################################################################################
    ####################################  -  SETUP DATA  -  ####################################
    
    params1, params2, params3, dataframes = read_data_init()

    nt, nx, ny, n, _ = params1
    T, dt, h, L, H, lbox, din, dbot = params2
    swing_start, pert_start, pert_dt, alpha, strouhal = params3
    
    cmap1, cmap2, cmap3, cmap4 = "Spectral_r", "bwr", "viridis", "turbo_r"
    strm_levels, strm_color, strm_alpha = 17, "grey", 0.5
    nStreakLines, nParticles, strk_lw, strk_color, strk_alpha = 4, 100, 3., "grey", 0.5

    
    # Oscillation parameters
    t_array = np.linspace(0., T, nt + 1)
    kappa = alpha / (2 * np.pi * strouhal)
    delta_x = 0. * kappa * (1. - np.cos(2. * np.pi * strouhal * t_array))
    delta_y = 0.5 * pert_dt / (2. * np.pi) * (1. - np.cos(2. * np.pi * (t_array - pert_start) / pert_dt))
    delta_x[t_array < swing_start] = 0.
    delta_y[~((pert_start < t_array) * (t_array < pert_start + pert_dt))] = 0.
    
    # Mesh points
    xm, ym  = np.linspace(h/2., L-h/2., nx), np.linspace(h/2., H-h/2., ny)  # middle of cell : p value of MAC mesh
    x, y = np.linspace(0, L, nx + 1), np.linspace(0, H, ny + 1)  # corners of cell : w value of MAC mesh
    (xx_, yy_), (xxm, yym) = np.meshgrid(x, y), np.meshgrid(xm, ym)
    mask_middle = (din < xxm) * (xxm < din+lbox) * (dbot < yym) * (yym < dbot+1)
    mask_corner = (din <= xx_) * (xx_ <= din+lbox) * (dbot <= yy_) * (yy_ <= dbot+1.)

    # Fields
    p, u, v, w, p_masked, w_masked, psi = read_block(dataframes)


    #############################################################################################
    ######################################  -  FIGURE  -  #######################################
    fig, axs = plt.subplots(2, 1, figsize=(12., 8.), constrained_layout=False, sharex="all", sharey="all")

    # pmin, pmax = np.amin(p), np.amax(p)
    pmin, pmax = -2.25, 2.25
    # wmin, wmax = np.amin(w)/10., np.amax(w)/10.
    wmin, wmax = -25., 25.

    # pressure and vorticity fields
    pressure = axs[0].imshow(np.zeros_like(p), extent=(0, L, 0, H), vmin=pmin, vmax=pmax, cmap=cmap1, origin="lower")
    vorticity = axs[1].imshow(np.zeros_like(w), extent=(0, L, 0, H), vmin=wmin, vmax=wmax, cmap=cmap2, origin="lower")
    cax_p, cax_w = make_colorbar_with_padding(axs[0]), make_colorbar_with_padding(axs[1])
    cbar_p, cbar_w = fig.colorbar(pressure, cax=cax_p), fig.colorbar(vorticity, cax=cax_w)

    # streamlines and streaklines
    if kwargs["stream"]:
        strm = [axs[0].contour(xx_, yy_, psi, colors=strm_color, alpha=strm_alpha, levels=strm_levels)]
    if kwargs["streak"]:
        spots_y = np.linspace(H/(nStreakLines+1), H - H/(nStreakLines+1), nStreakLines)
        streaklines_x = [-np.ones(nParticles) for _ in range(nStreakLines)]
        streaklines_y = [spots_y[j]*np.ones(nParticles) for j in range(nStreakLines)]
        lines = [axs[1].plot([], [], "-", marker=".", lw=strk_lw, color=strk_color, markersize=3,
                             alpha=strk_alpha)[0] for _ in range(nStreakLines)]

    # Obstacle
    patches = [Rectangle((din, dbot), lbox, 1., fc="grey", ec="none", alpha=1., zorder=5) for _ in axs]
    
    # Axis limits
    fig.tight_layout()
    for ax in axs:
        ax.axis([0., L, 0. + 0.*H, H - 0.*H])
        # ax.axis([2., 10., 1.5, 3.5])
        ax.set_aspect('equal')
    
    # Time text
    time_str = "t = {:5.3f} [s]"
    bbox_dic = dict(boxstyle="round", fc="wheat", ec="none", alpha=0.85)
    time_text = axs[0].text(0.8,0.9, time_str.format(0), fontsize=ftSz2, transform=axs[0].transAxes, bbox=bbox_dic)


    save = "none"
    # save = "gif"
    # save = "mp4"
    # save = "html"

    if save == "none":
        # init()
        # for i in range(nt//10):
        #     update(i)
        # plt.show()

        anim = FuncAnimation(fig, update, tqdm(range(nt//10+1)), interval=100, blit=False, init_func=init, repeat=False)
        plt.show()
    
    elif save == "gif":
        anim = FuncAnimation(fig, update, tqdm(range(nt+1)), interval=100, blit=False, init_func=init, repeat=False)
        writerGIF = PillowWriter(fps=15)
        anim.save("flow.gif", writer=writerGIF)
    
    elif save == "mp4":
        anim = FuncAnimation(fig, update, tqdm(range(nt+1)), interval=100, blit=False, init_func=init, repeat=False)
        writerMP4 = FFMpegWriter(fps=15)
        anim.save("flow.mp4", writer=writerMP4)

    elif save == "html":
        fig.subplots_adjust(bottom=0.02, top=0.98, left=0.02, right=0.98, hspace=0.05)
        init()
        for t in tqdm(range(nt + 1)):
            update(t)
            fig.savefig("./anim/frame_{:05d}.png".format(t), format="png", bbox_inches='tight', pad_inches=0.02)
        
        modify_html("animation.html")
