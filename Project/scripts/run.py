import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import sys

from tqdm import tqdm
from scipy.interpolate import RegularGridInterpolator
from time import perf_counter
from matplotlib.animation import FuncAnimation, PillowWriter, FFMpegWriter
from matplotlib.patches import Rectangle
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import checkdep_usetex


class Simulation:

    def __init__(self, filenames, analyze=False):

        params1, params2, params3, params4, params5, dataframes, diagnostics, self.u_avg, self.v_avg = read_data_init(filenames, analyze)

        self.nx, self.ny, self.n, self.temperature = params1
        self.RE, self.T, self.dt, self.h, self.save_freq, self.start_avg = params2
        self.L, self.H, self.lbox, self.din, self.dbot = params3
        self.alpha, self.st, self.swing_start, self.kappa_y, self.st_y, self.pert_start, self.smooth, self.n_cycles = params4
        self.prandtl, self.grashof, self.eckert, self.tmin, self.tmax = params5

        self.df_p, self.df_u, self.df_v = dataframes[:3]
        self.df_T = dataframes[-1]

        # self.nt_simu = self.nt
        # self.dt_simu = self.dt
        # self.dt *= self.save_modulo
        # self.nt //= self.save_modulo
        # self.t_simu = np.linspace(0., self.T, self.nt_simu + 1)
        # self.t = np.linspace(0., self.T, self.nt + 1)

        self.t_simu = diagnostics[:, 6]
        self.t = self.t_simu[diagnostics[:, 7] == 1]
        self.nt = self.t.shape[0] - 1

        self.kappa = self.alpha / (2. * np.pi * self.st)
        self.alpha_y = self.kappa_y * (2. * np.pi * self.st_y)
        self.pert_dt = self.n_cycles / self.st_y

        self.delta_x = self.kappa * (1. - np.cos(2. * np.pi * self.st * (self.t - self.swing_start)))
        self.u_mesh = self.alpha * np.sin(2. * np.pi * self.st * (self.t - self.swing_start))
        
        g1 = np.exp((self.pert_start - self.t) / self.smooth) if self.smooth > 1e-6 else 0.
        g2 = np.exp((self.t - self.pert_start - self.pert_dt) / self.smooth) if self.smooth > 1e-6 else 0.
        self.delta_y = self.kappa_y * (1. - np.cos(2. * np.pi * self.st_y * (self.t - self.pert_start))) * (1. - g1) * (1. - g2)
        if self.smooth > 1e-6:
            self.v_mesh = self.alpha_y * np.sin(2. * np.pi * self.st_y * (self.t - self.pert_start)) * g1 / self.smooth * (-g2) / self.smooth
        else:
            self.v_mesh = self.alpha_y * np.sin(2. * np.pi * self.st_y * (self.t - self.pert_start))

        self.delta_x[self.t < self.swing_start] = 0.
        self.delta_y[~((self.pert_start < self.t) * (self.t < self.pert_start + self.pert_dt))] = 0.
        self.u_mesh[self.t < self.swing_start] = 0.
        self.v_mesh[~((self.pert_start < self.t) * (self.t < self.pert_start + self.pert_dt))] = 0.

        self.xm = np.linspace(self.h/2., self.L-self.h/2., self.nx)   # middle of cell : p value of MAC mesh
        self.ym = np.linspace(self.h/2., self.H-self.h/2., self.ny)   
        self.x = np.linspace(0, self.L, self.nx + 1)  # corners of cell : w value of MAC mesh
        self.y = np.linspace(0, self.H, self.ny + 1)
        (self.xx, self.yy), (self.xxm, self.yym) = np.meshgrid(self.x, self.y), np.meshgrid(self.xm, self.ym)
        self.mask_middle = (self.din < self.xxm) * (self.xxm < self.din+self.lbox) * (self.dbot < self.yym) * (self.yym < self.dbot+1)
        self.mask_corner = (self.din <= self.xx) * (self.xx <= self.din+self.lbox) * (self.dbot <= self.yy) * (self.yy <= self.dbot+1)

        self.reh = diagnostics[:, 0]
        self.rew = diagnostics[:, 1]
        self.drag_p = diagnostics[:, 2]
        self.drag_f = diagnostics[:, 3]
        self.lift_p = diagnostics[:, 4]
        self.lift_f = diagnostics[:, 5]


def read_data_init(filenames, analyze=False):
    
    filename_params, filename_p, filename_u, filename_v, filename_T, filename_stats, filename_u_avg, filename_v_avg = filenames

    with open(filename_params, "r") as f:
        lines = f.readlines()
        params1 = [int(x) for x in lines[0].split(" ")]
        params2 = [float(x) for x in lines[1].split(" ")]
        params3 = [float(x) for x in lines[2].split(" ")]
        params4 = [float(x) for x in lines[3].split(" ")]
        params5 = [float(x) for x in lines[4].split(" ")]

    nx, ny, n, temp_enabled = params1
    for i in range(2, len(params3)):
        params3[i] = int(params3[i])

    size_p, size_u, size_v = nx * ny, (nx+1) * (ny+2), (nx+2) * (ny+1)
    
    df_p = pd.read_csv(filename_p, chunksize=size_p, header=None)
    df_u = pd.read_csv(filename_u, chunksize=size_u, header=None)
    df_v = pd.read_csv(filename_v, chunksize=size_v, header=None)
    dfs = [df_p, df_u, df_v]

    # if bool(params4[0]):
    if temp_enabled == 1:
        dfs.append(pd.read_csv(filename_T, chunksize=size_p, header=None))

    diagnostics = np.loadtxt(filename_stats)
    u_avg = np.array(pd.read_csv(filename_u_avg, header=None)).reshape(nx+1, ny+2)
    v_avg = np.array(pd.read_csv(filename_v_avg, header=None)).reshape(nx+2, ny+1)
    
    return params1, params2, params3, params4, params5, dfs, diagnostics, u_avg, v_avg
    
    # find bounds
    # t1 = perf_counter()
    # (pmin, pmax), (umin, umax), (vmin, vmax), (wmin, wmax) = find_all_bounds()
    # print(f"Time to find bound = {perf_counter() - t1:.3f} s")


# wb = [0., 0.]
def read_block(sim, needed={}):
    df_p, df_u, df_v = sim.df_p, sim.df_u, sim.df_v
    ret_list = []
    
    if needed.get("p", False):
        p = np.array(next(df_p)).reshape((sim.nx, sim.ny))
        ret_list.append(p)
    if needed.get("u", False) or needed.get("w", False):
        u = np.array(next(df_u)).reshape((sim.nx + 1, sim.ny + 2))
        if needed.get("u", False):
            ret_list.append(u)
    if needed.get("v", False) or needed.get("w", False):
        v = np.array(next(df_v)).reshape((sim.nx + 2, sim.ny + 1))
        if needed.get("v", False):
            ret_list.append(v)
    if needed.get("w", False):
        w = ((v[1:, :] - v[:-1, :]) - (u[:, 1:] - u[:, :-1])) / sim.h  # dv/dx - du/dy --> values at corners
        ret_list.append(w)
    if needed.get("T", False):
        T = np.array(next(sim.df_T)).reshape((sim.nx, sim.ny)) if sim.temperature else None
        ret_list.append(T)

    # wb[0] = min(np.amin(w), wb[0])
    # wb[1] = max(np.amax(w), wb[1])

    return ret_list


def apply_mask(sim, p, w, T):
    
    w_bis = 0.25 * (w[1:, 1:] + w[1:, :-1] + w[:-1, 1:] + w[:-1, :-1])  # average w on corners to get value in the middle
    w_masked = np.ma.masked_array(w_bis.T, sim.mask_middle)
    T_masked = np.ma.masked_array(T.T, sim.mask_middle) if sim.temperature else None
    p_masked = np.ma.masked_array(p.T, sim.mask_middle)

    return p_masked, w_masked, T_masked


def compute_psi(sim, u):
    # d(psi) = (+v) dx + (-u) dy

    psi = np.zeros((sim.nx+1, sim.ny+1))
    
    n = sim.n

    iwl, iwr = sim.din*n + 1, (sim.din + sim.lbox)*n + 1
    jwb, jwa = sim.dbot*n + 1, (sim.dbot + 1)*n + 1

    psi[:iwl, 1:] = sim.h * np.cumsum(-u[:iwl, 1:-1], axis=1)                 # block 1
    psi[iwl:iwr, 1:jwb] = sim.h * np.cumsum(-u[iwl:iwr, 1:jwb], axis=1)       # block 2
    psi[iwl:iwr, jwa-1] = psi[iwl:iwr, dbot*n] + np.sum(-u[iwl-1, jwb:jwa]) * sim.h       # block 3  # (-u_mesh[t]) * 1. 
    psi[iwl:iwr, jwa:] = psi[iwl:iwr, jwa-1, np.newaxis] + sim.h * np.cumsum(-u[iwl:iwr, jwa:-1], axis=1)
    psi[iwr:, 1:] = sim.h * np.cumsum(-u[iwr:, 1:-1], axis=1)  # block 4
    
    return psi


def find_bounds(filename):
    # min_value, max_value = np.inf, -np.inf
    min_value, max_value = 0., 0.
    with open(filename, "r") as f:
        for line in f:
            value = float(line)
            # min_value = min(min_value, value)
            # max_value = max(max_value, value)
            min_value = min(min_value, max(value, -100.))
            max_value = max(max_value, min(value, 100.))
    return min_value, max_value


def find_all_bounds(filename_p, filename_u, filename_v):
    res = [(0., 0.) for _ in range(3)]
    if filename_p != "":
        print("Computing bounds for p", end="", flush=True)
        bounds_p = find_bounds(filename_p)
        res[0] = bounds_p

    if filename_u != "":
        print("\rComputing bounds for u", end="")
        bounds_u = find_bounds(filename_u)
        res[1] = bounds_u

    if filename_v != "":
        print("\rComputing bounds for v", end="")
        bounds_v = find_bounds(filename_v)
        res[2] = bounds_v
    
    print("\rComputing bounds finished")
    return res


def update(t_idx):

    p, u, v, w, T = read_block(sim, needed={"p": True, "u": True, "v": True, "w": True, "T":True})
    p_masked, w_masked, T_masked = apply_mask(sim, p, w, T)

    # print(p[0, 0], p[0, 1], p[1, 0])

    animated_items = []

    for patch in patches:
        patch.set_xy([sim.din + sim.delta_x[t_idx], sim.dbot + sim.delta_y[t_idx]])
        animated_items.append(patch)
     
    pressure.set_extent([sim.delta_x[t_idx], sim.L + sim.delta_x[t_idx], 0 + sim.delta_y[t_idx], sim.H + sim.delta_y[t_idx]])  # update position
    pressure.set_data(p_masked)
    animated_items.append(pressure)

    vorticity.set_extent([sim.delta_x[t_idx], sim.L + sim.delta_x[t_idx], 0 + sim.delta_y[t_idx], sim.H + sim.delta_y[t_idx]])  # update position
    vorticity.set_data(w_masked)
    animated_items.append(vorticity)

    if sim.temperature:
        temperature.set_extent([sim.delta_x[t_idx], sim.L + sim.delta_x[t_idx], 0 + sim.delta_y[t_idx], sim.H + sim.delta_y[t_idx]])  # update position
        temperature.set_data(T_masked)
        animated_items.append(temperature)


    if kwargs["stream"] and t_idx > 0:
        for item in strm[0].collections:
            item.remove()
        
        psi = compute_psi(sim, u)
        psi = np.ma.masked_array(psi.T, sim.mask_corner)

        psi_min, psi_max = np.amin(psi), np.amax(psi)
        psi_avg, psi_dif = (psi_max + psi_min) / 2., psi_max - psi_min
        levels = np.linspace(psi_min + 0.01 * psi_dif, psi_max - 0.01 * psi_dif, nStreamLines)
        levels = levels - scaling * psi_dif / (2. * np.pi) * np.sin(2. * np.pi / psi_dif * (levels - psi_avg))

        strm[0] = axs[0].contour(sim.xx + sim.delta_x[t_idx], sim.yy + sim.delta_y[t_idx], psi, linewidths=strm_lw,
                                 colors=strm_color, alpha=strm_alpha, levels=levels, linestyles="solid")
        animated_items.append(strm[0].collections)
    
    if kwargs["streak"]:
        U, V = 0.5 * (u[:, :-1] + u[:, 1:]), 0.5 * (v[:-1, :] + v[1:, :])
        u_interp = RegularGridInterpolator((sim.x, sim.y), U, bounds_error=False, fill_value=0)
        v_interp = RegularGridInterpolator((sim.x, sim.y), V, bounds_error=False, fill_value=0)
        cycle = int(np.ceil((L / 1.) / nParticles / sim.dt)) # [temps pour traverser / nbre particules] = delai
        
        for j, (line, x_particles, y_particles) in enumerate(zip(lines, streaklines_x, streaklines_y)):
            # if t % cycle == 0:
            x_particles[(t_idx//cycle) % nParticles] = 2 * sim.kappa
            y_particles[(t_idx//cycle) % nParticles] = spots_y[j]
            mask = (0. <= x_particles) * (x_particles <= L) * (0. <= y_particles) * (y_particles <= H)
            
            # evaluate the field at the right place: take the shift in consideration
            args = np.c_[x_particles[mask] - sim.delta_x[t_idx], y_particles[mask] - sim.delta_y[t_idx]]
            
            x_particles[mask] += sim.dt * u_interp(args)  # Euler explicite
            y_particles[mask] += sim.dt * v_interp(args)  # Euler explicite

            x_particles[~mask] = np.nan
            y_particles[~mask] = np.nan
            line.set_data(x_particles, y_particles)
        
        animated_items += lines


    time_text.set_text(time_str.format(sim.t[t_idx]))
    
    idx_t_simu = np.argwhere(sim.t_simu == sim.t[t_idx])[0][0]
    idx_t_simu = idx_t_simu if t_idx < nt else idx_t_simu - 1
    dt_now = sim.t_simu[idx_t_simu+1] - sim.t_simu[idx_t_simu]
    mantisse_dt = dt_now / np.power(10., -3)
    # dt_text = r"    $\Delta t^*= {:.1f}\cdot 10^{{{:.0f}}}$".format(mantisse_dt, -3)
    axs[0].set_title(title_text + dt_text.format(mantisse_dt, -3), fontsize=title_fontsize)

    if save != "html":
        pbar.update(1)

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


def create_html():
    html_filename = f"{path_anim:s}.html"
    line_nb = 152
    with open(html_ref, "r") as f:
        html_lines = f.readlines()
        new_line = f"  call('{case_name:s}', {sim.nt+1:d})\n"
        html_lines[line_nb] = new_line

    with open(html_filename, "w") as f:
        f.write("".join(html_lines))
        f.close()


ftSz1, ftSz2, ftSz3 = 20, 17, 15  # time text, parameters, labels
plt.rcParams['font.family'] = 'monospace'
plt.rcParams["text.usetex"] = checkdep_usetex(True)
fix_title = True

if __name__ == "__main__":

    kwargs = {"stream":True, "streak":False}
    
    # case_name = "no_slip_gr"  # case_1
    if (len(sys.argv) != 5) or ("-h" in sys.argv):
        print("usage: python3 run.py -dir [directory_name] -save [no, gif, mp4, html]") 
    if sys.argv[1] == "-dir":
        case_name = sys.argv[2]
    if (sys.argv[3] == "-save") and (sys.argv[4] in ["none", "gif", "mp4", "html"]):
        save = sys.argv[4]

    # path_src = "/media/vincelinux/LaCie/LMECA2660"
    path_src = ".."
    
    # path_dst = "/media/vincelinux/LaCie/LMECA2660"
    path_dst = ".."

    path_res = path_src + "/results/" + case_name
    path_anim = path_dst + "/anim/" + case_name
    html_ref = path_dst + "/anim/no_slip.html"  # copy and create new html with a different path to images

    filename_params = f"{path_res}/simu_params.txt"
    filename_p = f"{path_res}/simu_p.txt"
    filename_u = f"{path_res}/simu_u.txt"
    filename_v = f"{path_res}/simu_v.txt"
    filename_T = f"{path_res}/simu_T.txt"
    filename_stats = f"{path_res}/simu_stats.txt"
    filename_u_avg = f"{path_res}/simu_u_avg.txt"
    filename_v_avg = f"{path_res}/simu_v_avg.txt"

    sim = Simulation([filename_params, filename_p, filename_u, filename_v, filename_T,
                      filename_stats, filename_u_avg, filename_v_avg], analyze=False)

    cmap1, cmap2, cmap3, cmap4 = "Spectral_r", "bwr", "RdBu_r", "turbo_r"  # OrRd
    nStreamLines, scaling, strm_lw, strm_color, strm_alpha = 24, 0.70, 1.25, "grey", 0.25
    nStreakLines, nParticles, strk_lw, strk_color, strk_alpha = 2, 500, 3., "grey", 0.5


    # make some parameters global for access simplicity

    nt, nx, ny = sim.nt, sim.nx, sim.ny
    L, H, din, dbot, lbox = sim.L, sim.H, sim.din, sim.dbot, sim.lbox
    xx, yy = sim.xx, sim.yy


    #############################################################################################
    ######################################  -  FIGURE  -  #######################################

    n_plots = 3 if sim.temperature else 2
    height = 9. if sim.temperature else 8.
    width = (height * 0.9 / n_plots) * (sim.L + 2) / sim.H
    # width, height = (9., 8.5) if sim.temperature else (12., 8.)

    fig, axs = plt.subplots(n_plots, 1, figsize=(width, height), constrained_layout=False, sharex="all", sharey="all")

    # pmin, pmax = np.amin(p), np.amax(p)
    pmin, pmax = -2.75, 2.75
    # wmin, wmax = np.amin(w)/10., np.amax(w)/10.
    wmin, wmax = -22., 22.
    # Tmin, Tmax = find_bounds(filename_T)
    Tmin, Tmax = -1., 1.
    # Tmin, Tmax = 0., 0.016

    # pressure and vorticity fields
    pressure = axs[0].imshow(np.zeros((nx, ny)), extent=(0, L, 0, H), vmin=pmin, vmax=pmax, cmap=cmap1, origin="lower")
    vorticity = axs[1].imshow(np.zeros((nx+1, ny+1)), extent=(0, L, 0, H), vmin=wmin, vmax=wmax, cmap=cmap2, origin="lower")
    cax_p, cax_w = make_colorbar_with_padding(axs[0]), make_colorbar_with_padding(axs[1])
    cbar_p, cbar_w = fig.colorbar(pressure, cax=cax_p), fig.colorbar(vorticity, cax=cax_w)
    cax_p.set_ylabel(r"$\Delta p \: / \: \rho U^2$", fontsize=ftSz3)
    cax_w.set_ylabel(r"$\omega \: H_{{box}} \: / \: U$", fontsize=ftSz3)

    space_sep = r"\qquad\quad" if sim.eckert == 0. else r"\qquad"
    title_text = r"$Re={:.0f} {:s}$".format(sim.RE, space_sep)
    title_fontsize = ftSz2 * 15./17. if sim.temperature else ftSz2

    if sim.temperature:
        temperature = axs[2].imshow(np.zeros((nx, ny)), extent=(0, L, 0, H), vmin=Tmin, vmax=Tmax, cmap=cmap3, origin="lower")
        cax_T = make_colorbar_with_padding(axs[2])
        cbar_T = fig.colorbar(temperature, cax=cax_T)
        cax_T.set_yticks([Tmin + i * (Tmax - Tmin) / (5. - 1.) for i in range(5)])
        cax_T.set_ylabel(r"$(T - T_0)\: / \: \Delta T$", fontsize=ftSz3)

        title_text += "$Pr={:.1f} {:s}$".format(sim.prandtl, space_sep)
        if sim.grashof > 0:
            base10power_Gr = np.floor(np.log10(sim.grashof))
            mantisse_Gr = sim.grashof / np.power(10, base10power_Gr)
            title_text += "$Gr={:.0f} \cdot 10^{{{:.0f}}} {:s}$".format(mantisse_Gr, base10power_Gr, space_sep)
        else:
            title_text += "$Gr=0 {:s}$".format(space_sep)

        if sim.eckert > 0:
            base10power_Ec = np.floor(np.log10(sim.eckert))
            mantisse_Ec = sim.eckert / np.power(10, base10power_Ec)
            title_text += "$Ec={:.0f} \cdot 10^{{{:.0f}}} {:s}$".format(mantisse_Ec, base10power_Ec, space_sep)

    title_text += "$h^* = {:.3f} {:s}$".format(sim.h, space_sep)
    base10power_dt = np.floor(np.log10(sim.dt))
    mantisse_dt = sim.dt / np.power(10, base10power_dt)
    if fix_title:  # avoid the title to move when dt changes
        dt_text = r"$\Delta t^*= \mathtt{{{:.1f}}} \cdot 10^{{{:.0f}}}$"
    else:
        dt_text = r"$\Delta t^*= {:.1f}\cdot 10^{{{:.0f}}}$"

    axs[0].set_title(title_text + dt_text.format(mantisse_dt, base10power_dt), fontsize=title_fontsize)

    # streamlines and streaklines
    if kwargs["stream"]:
        strm = [axs[0].contour(xx, yy, xx*yy, colors=strm_color, alpha=0., levels=1)] # on purpose
    if kwargs["streak"]:
        spots_y = np.linspace(H/(nStreakLines+1), H - H/(nStreakLines+1), nStreakLines)
        streaklines_x = [-np.ones(nParticles) for _ in range(nStreakLines)]
        streaklines_y = [spots_y[j]*np.ones(nParticles) for j in range(nStreakLines)]
        lines = [axs[1].plot([], [], ls="", marker=".", lw=strk_lw, color=strk_color, markersize=3,
                             alpha=strk_alpha)[0] for _ in range(nStreakLines)]

    # Obstacle
    patches = [Rectangle((din, dbot), lbox, 1., fc="lightgrey", ec="none", alpha=1., zorder=5) for _ in axs]
    for ax, patch in zip(axs, patches):
        ax.add_patch(patch)
    
    # Axis limits
    fig.tight_layout()
    for ax in axs:
        ax.axis([0., L, 0., H])
        ax.set_aspect('equal')

    # Time text
    position_text = (0.03,0.85) if n_plots == 3 else (0.855, 0.9)
    fontsize = ftSz1 if plt.rcParams["text.usetex"] else ftSz2
    time_str = r"$t^*$ = {:5.2f}"
    bbox_dic = dict(boxstyle="round", fc="wheat", ec="none", alpha=0.85)
    time_text = axs[0].text(*position_text, time_str.format(0), fontsize=fontsize,
                             transform=axs[0].transAxes, bbox=bbox_dic)


    ##############################################################################################
    ###################################  -  SAVE & DISPLAY  -  ###################################

    # save = "none"
    # save = "gif"
    # save = "mp4"
    # save = "html"

    if save == "none":
        pbar = tqdm(total=nt+1)
        anim = FuncAnimation(fig, update, nt+1, interval=40, blit=False, init_func=lambda :None, repeat=False)
        pbar.close()
        plt.show()
    
    elif save == "gif":
        pbar = tqdm(total=nt+1)
        anim = FuncAnimation(fig, update, nt+1, interval=40, blit=False, init_func=lambda :None, repeat=False)
        writerGIF = PillowWriter(fps=15)
        anim.save(f"{path_anim}.gif", writer=writerGIF, dpi=100)
        pbar.close()
    
    elif save == "mp4":
        # fig.subplots_adjust(left=0.0)
        pbar = tqdm(total=nt+1)
        anim = FuncAnimation(fig, update, nt+1, interval=40, blit=False, init_func=lambda :None, repeat=False)
        writerMP4 = FFMpegWriter(fps=25)
        anim.save(f"{path_anim}.mp4", writer=writerMP4)
        pbar.close()

    elif save == "html":
        os.makedirs(path_anim + "/", exist_ok=True)
        fig.subplots_adjust(bottom=0.01, top=0.99, left=0.01, right=0.99, hspace=0.05)
        for t in tqdm(range(nt+1)):
            update(t)
            fig.savefig(f"{path_anim:s}/frame_{t:05d}.png", format="png", bbox_inches='tight', pad_inches=0.055)
        
        create_html()
