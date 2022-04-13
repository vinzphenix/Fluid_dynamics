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

path_dir = "../results"
filename_params = f"{path_dir}/simu_params.txt"
filename_p = f"{path_dir}/simu_p.txt"
filename_u = f"{path_dir}/simu_u.txt"
filename_v = f"{path_dir}/simu_v.txt"
global size_p, size_u, size_v

def read_data_init():

    with open(filename_params, "r") as f:
        lines = f.readlines()
        params1 = [int(x) for x in lines[0].split(" ")]
        params2 = [float(x) for x in lines[1].split(" ")]
        params3 = [float(x) for x in lines[2].split(" ")]
    
    nt, nx, ny, n, save_modulo = params1
    params2[2] *= save_modulo  # dt multiplied since not all times are saved

    for i in range(6, 9):
        params2[i] = int(params2[i])

    global size_p, size_u, size_v
    size_p, size_u, size_v = nx * ny, (nx+1) * (ny+2), (nx+2) * (ny+1)
    
    df_p = pd.read_csv(filename_p, chunksize=size_p, header=None)
    df_u = pd.read_csv(filename_u, chunksize=size_u, header=None)
    df_v = pd.read_csv(filename_v, chunksize=size_v, header=None)
    
    # find bounds
    # t1 = perf_counter()
    # (pmin, pmax), (umin, umax), (vmin, vmax), (wmin, wmax) = find_all_bounds()
    # print(f"Time to find bound = {perf_counter() - t1:.3f} s")

    return params1, params2, params3, [df_p, df_u, df_v]


def compute_psi(u, v, t):
    # d(psi) = (+v) dx + (-u) dy

    psi = np.zeros((nx+1, ny+1))
    # uu = u - u_mesh[t]
    # vv = v - v_mesh[t]

    iwl, iwr = din*n+1, (din+lbox)*n+1
    jwb, jwa = dbot*n+1, (dbot+1)*n+1

    psi[:iwl, 1:] = h * np.cumsum(-u[:iwl, 1:-1], axis=1)                 # block 1
    psi[iwl:iwr, 1:jwb] = h * np.cumsum(-u[iwl:iwr, 1:jwb], axis=1)       # block 2
    psi[iwl:iwr, jwa-1] = psi[iwl:iwr, dbot*n] + np.sum(-u[iwl-1, jwb:jwa]) * h       # block 3  # (-u_mesh[t]) * 1. 
    psi[iwl:iwr, jwa:] = psi[iwl:iwr, jwa-1, np.newaxis] + h * np.cumsum(-u[iwl:iwr, jwa:-1], axis=1)
    psi[iwr:, 1:] = h * np.cumsum(-u[iwr:, 1:-1], axis=1)  # block 4
    
    return psi


def read_block(dataframe_list):
    df_p, df_u, df_v = dataframe_list
    
    p = np.array(next(df_p)).reshape((nx, ny))
    u = np.array(next(df_u)).reshape((nx+1, ny+2))
    v = np.array(next(df_v)).reshape((nx+2, ny+1))
    
    w = ((v[1:, :] - v[:-1, :]) - (u[:, 1:] - u[:, :-1])) / h  # dv/dx - du/dy --> values at corners

    u = (u[:, 1:] + u[:, :-1]) / 2.
    v = (v[1:, :] + v[:-1, :]) / 2.

    return p, u, v, w


def make_colorbar_with_padding(ax):
    """
    Create colorbar axis that fits the size of a plot - detailed here: http://chris35wills.github.io/matplotlib_axis/
    """
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="2%", pad=0.15)
    return(cax) 


def plot_max_RE():
    reh = np.empty(nt+1)
    rew = np.empty(nt+1)
    t_array = np.linspace(0, T, nt+1)

    for t in range(nt+1):
        p, u, v, w = read_block(dataframes)
        reh[t] = np.amax(np.abs(u[~mask_corner]) + np.abs(v[~mask_corner])) * h * RE
        rew[t] = np.amax(np.abs(w[~mask_corner])) * h * h * RE

    fig, ax = plt.subplots(1, 1, figsize=(10., 6.))
    fig.tight_layout()
    ax.plot(t_array, reh, marker=".", markersize=5, color='C0', label='RE_h')
    ax.plot(t_array, rew, marker=".", markersize=5, color='C1', label='RE_w')

    ax.set_xlabel(r"$t U_{\infty} / H_{box}$", fontsize=ftSz2)
    ax.set_ylabel(r"$Re$", fontsize=ftSz2)
    ax.legend(fontsize=ftSz3)
    ax.grid(ls=':')
    plt.show()

    return


def plot_lift_drag():
    drag = np.zeros(nt+1)
    lift = np.zeros(nt+1)

    return


def plot_streamlines(t_list):
    fig, axs = plt.subplots(2, 1, figsize=(10., 6.), sharex="all", sharey="all")
    fig.tight_layout()

    i = 0
    for t in range(nt+1):
        if t not in t_list:
            continue

        p, u, v, w = read_block(dataframes)
        axs[i].streamplot(x + delta_x[t], y + delta_y[t], u.T, v.T, density=1., color=np.hypot(u.T, v.T), cmap=cmap)
        i += 1

    bbox_dic = dict(boxstyle="round", fc="wheat", ec="none", alpha=0.85)

    axs[-1].set_xlabel(r"$x / H_{box}$", fontsize=ftSz2)
    for t, ax in zip(t_list, axs):
        box = Rectangle((din, dbot), lbox, 1., fc="grey", ec="none", alpha=1., zorder=5)
        ax.add_patch(box)

        ax.text(0.8,0.85, "t = {:5.2f} s".format(t * dt), fontsize=ftSz2, transform=ax.transAxes, bbox=bbox_dic)
        ax.set_ylabel(r"$y / H_{{box}}$", fontsize=ftSz2)
        ax.axis([0., L, 0., H])

    plt.show()
    return


def plot_vorticity(t_list):
    fig, axs = plt.subplots(2, 1, figsize=(10., 6.), sharex="all", sharey="all")
    fig.tight_layout()
    caxs = make_colorbar_with_padding(axs[0]), make_colorbar_with_padding(axs[1])
    contours = [None, None]
    levels = np.linspace(-25., 25., 10)
    
    i = 0
    for t in range(nt+1):
        if t not in t_list:
            continue

        p, u, v, w = read_block(dataframes)
        contours[i] = axs[i].contour(x + delta_x[t], y + delta_y[t], w.T, levels, cmap=cmap)
        cbar = fig.colorbar(contours[i], cax=caxs[i])
        i += 1

    bbox_dic = dict(boxstyle="round", fc="wheat", ec="none", alpha=0.85)

    axs[-1].set_xlabel(r"$x / H_{box}$", fontsize=ftSz2)
    for t, ax in zip(t_list, axs):
        box = Rectangle((din + delta_x[t], dbot + delta_y[t]), lbox, 1., fc="grey", ec="none", alpha=1., zorder=5)
        ax.add_patch(box)

        ax.text(0.8,0.85, "t = {:5.2f} s".format(t * dt), fontsize=ftSz2, transform=ax.transAxes, bbox=bbox_dic)
        ax.set_ylabel(r"$y / H_{{box}}$", fontsize=ftSz2)
        ax.axis([0., L, 0., H])

    plt.show()
    return


def plot_average_flow():
    fig, axs = plt.subplots(2, 1, figsize=(10., 6.), sharex="all", sharey="all")
    fig.tight_layout()
    caxs = make_colorbar_with_padding(axs[0]), make_colorbar_with_padding(axs[1])
    
    n_start = 3

    p_sum = np.zeros((nx, ny))
    u_sum = np.zeros((nx+1, ny+1))
    v_sum = np.zeros((nx+1, ny+1))
    for t in range(n_start, nt+1):

        p, u, v, w = read_block(dataframes)
        p_sum += p
        u_sum += u
        v_sum += v

    p_sum /= (nt+1 - n_start)
    u_sum /= (nt+1 - n_start)
    v_sum /= (nt+1 - n_start)

    psi = axs[0].contour(x + delta_x[t], y + delta_y[t], u_sum.T, 10, cmap=cmap)
    w = axs[1].contour(x + delta_x[t], y + delta_y[t], v_sum.T, 10, cmap=cmap)
    cbar1 = fig.colorbar(psi, cax=caxs[0])
    cbar2 = fig.colorbar(w, cax=caxs[1])

    bbox_dic = dict(boxstyle="round", fc="wheat", ec="none", alpha=0.85)

    axs[0].set_title(r"$\psi$", fontsize=ftSz2)
    axs[1].set_title(r"$\omega$", fontsize=ftSz2)
    axs[-1].set_xlabel(r"$x / H_{box}$", fontsize=ftSz2)
    for ax in axs:
        ax.text(0.8,0.85, "t = {:5.2f} s".format(t * dt), fontsize=ftSz2, transform=ax.transAxes, bbox=bbox_dic)
        ax.set_ylabel(r"$y / H_{{box}}$", fontsize=ftSz2)
        ax.axis([0., L, 0., H])

    plt.show()
    return



if __name__ == "__main__":

    cmap = plt.get_cmap("viridis")
    params1, params2, params3, dataframes = read_data_init()

    nt, nx, ny, n, _ = params1
    RE, T, dt, h, L, H, lbox, din, dbot = params2
    swing_start, pert_start, pert_dt, alpha, strouhal = params3

    t_array = np.linspace(0., T, nt + 1)
    kappa = alpha / (2 * np.pi * strouhal)
    delta_x = kappa * (1. - np.cos(2. * np.pi * strouhal * (t_array - swing_start)))
    delta_y = 0.5 * pert_dt / (2. * np.pi) * (1. - np.cos(2. * np.pi * (t_array - pert_start) / pert_dt))    
    u_mesh = alpha * np.sin(2. * np.pi * strouhal * (t_array - swing_start))
    v_mesh = 0.5 * np.sin(2. * np.pi * (t_array - pert_start) / pert_dt)
    delta_x[t_array < swing_start] = 0.
    delta_y[~((pert_start < t_array) * (t_array < pert_start + pert_dt))] = 0.
    u_mesh[t_array < swing_start] = 0.
    v_mesh[~((pert_start < t_array) * (t_array < pert_start + pert_dt))] = 0.

    delta_x += 1

    xm, ym  = np.linspace(h/2., L-h/2., nx), np.linspace(h/2., H-h/2., ny)  # middle of cell : p value of MAC mesh
    x, y = np.linspace(0, L, nx + 1), np.linspace(0, H, ny + 1)  # corners of cell : w value of MAC mesh
    (xx_, yy_), (xxm, yym) = np.meshgrid(x, y, indexing="ij"), np.meshgrid(xm, ym, indexing="ij")
    mask_middle = (din < xxm) * (xxm < din+lbox) * (dbot < yym) * (yym < dbot+1)
    mask_corner = (din < xx_) * (xx_ < din+lbox) * (dbot < yy_) * (yy_ < dbot+1.)

    # plot_max_RE()
    # plot_lift_drag()
    # plot_streamlines([1, 3])
    # plot_vorticity([5, 8])
    plot_average_flow()
