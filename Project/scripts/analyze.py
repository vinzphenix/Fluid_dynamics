# import os
# sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), 'run.py')))

from run import *


def plot_vorticity(t_list, cmap):

    n_plots = len(t_list)
    indices = [np.argmin(np.abs(this_t - sim.t)) for this_t in t_list]
        
    fig, axs = plt.subplots(n_plots, 1, figsize=(10., 10.), sharex="all", sharey="all")
    fig.tight_layout()
    caxs = [make_colorbar_with_padding(axs[i]) for i in range(n_plots)]
    
    i = 0
    for t in tqdm(range(nt+1)):
        p, u, v, w, _, w_masked = read_block(sim)
        if t not in indices:
            continue

        w = np.ma.masked_array(w.T, mask_corner)
        
        wmin, wmax = np.amin(w[~mask_hard]), np.amax(w[~mask_hard])
        wmax = max(np.abs(wmin), np.abs(wmax))
        # levels = np.linspace(wmax * 0.05, np.round(wmax * 0.95, 0), 10)
        # levels = np.r_[-levels[::-1], levels]
        levels = np.linspace(wmin, wmax, 20)
        # levels = levels - 0.5 * (2*wmax) / (2. * np.pi) * np.sin(2. * np.pi / (2 * wmax) * (levels))

        ct = axs[i].contour(x + delta_x[t], y + delta_y[t], w, levels=levels, cmap=cmap)
        cbar = fig.colorbar(ct, cax=caxs[i], ticks=np.round(np.linspace(levels[0], levels[-1], 7)))
        i += 1

    bbox_dic = dict(boxstyle="round", fc="wheat", ec="none", alpha=0.85)
    axs[-1].set_xlabel(r"$x / H_{box}$", fontsize=ftSz2)
    for idx, ax in zip(indices, axs):
        box = Rectangle((din + delta_x[idx], dbot + delta_y[idx]), lbox, 1., fc="grey", ec="none", alpha=1., zorder=5)
        ax.add_patch(box)
        ax.text(0.8,0.85, "t = {:5.2f} s".format(idx * sim.dt), fontsize=ftSz2, transform=ax.transAxes, bbox=bbox_dic)
        ax.set_ylabel(r"$y / H_{{box}}$", fontsize=ftSz2)
        ax.axis([0., L, 0., H])
        ax.set_aspect("equal")

    plt.show()
    return


def plot_streamlines(t_list, cmap):

    n_plots = len(t_list)
    indices = [np.argmin(np.abs(this_t - sim.t)) for this_t in t_list]
    x, y, delta_x, delta_y = sim.x, sim.y, sim.delta_x, sim.delta_y

    fig, axs = plt.subplots(n_plots, 1, figsize=(10., 10.), sharex="all", sharey="all")
    fig.tight_layout()

    i = 0
    for idx in tqdm(range(nt+1)):
        p, u, v, w, _, _ = read_block(sim)
        if idx not in indices:
            continue
        u = (u[:, 1:] + u[:, :-1]) / 2.
        v = (v[1:, :] + v[:-1, :]) / 2.
        speed = np.hypot(u.T, v.T)
        axs[i].streamplot(x + delta_x[idx], y + delta_y[idx], u.T, v.T, density=1.25, color=speed, cmap=cmap)
        i += 1

    bbox_dic = dict(boxstyle="round", fc="wheat", ec="none", alpha=0.85)
    axs[-1].set_xlabel(r"$x / H_{box}$", fontsize=ftSz2)
    for idx, ax in zip(indices, axs):
        box = Rectangle((din + delta_x[idx], dbot + delta_y[idx]), lbox, 1., fc="lightgrey", ec="none", alpha=1., zorder=5)
        ax.add_patch(box)
        ax.text(0.8,0.85, "t = {:5.2f} s".format(idx * sim.dt), fontsize=ftSz2, transform=ax.transAxes, bbox=bbox_dic)
        ax.set_ylabel(r"$y / H_{{box}}$", fontsize=ftSz2)
        ax.axis([0., L, 0., H])
        ax.set_aspect("equal")

    plt.show()
    return


def plot_average_flow(t_start, cmap1, cmap2):
    fig, axs = plt.subplots(2, 1, figsize=(10., 8.), sharex="all", sharey="all")
    fig.tight_layout()
    caxs = [make_colorbar_with_padding(axs[i]) for i in range(2)]

    idx_start = np.argmin(np.abs(sim.t - t_start))

    p_avg = np.zeros((nx, ny))
    u_avg = np.zeros((nx+1, ny+2))
    v_avg = np.zeros((nx+2, ny+1))

    for idx in tqdm(range(nt + 1)):
        p, u, v, w, p_masked, w_masked = read_block(sim)
        
        if idx < idx_start:
            continue
        
        p_avg += p
        u_avg += u
        v_avg += v

    p_avg /= (nt+1 - idx_start)
    u_avg /= (nt+1 - idx_start)
    v_avg /= (nt+1 - idx_start)

    w = ((v_avg[1:, :] - v_avg[:-1, :]) - (u_avg[:, 1:] - u_avg[:, :-1])) / sim.h
    w = np.ma.masked_array(w.T, mask_corner)

    u_avg = (u_avg[:, 1:] + u_avg[:, :-1]) / 2.
    v_avg = (v_avg[1:, :] + v_avg[:-1, :]) / 2.

    # STREAMLINES
    # psi = axs[0].contour(x, y, u_sum.T, 10, cmap=cmap)
    # cbar1 = fig.colorbar(psi, cax=caxs[0])
    speed = np.hypot(u_avg.T, v_avg.T)
    axs[0].streamplot(x, y, u_avg.T, v_avg.T, density=1.25, color=speed, cmap=cmap1)
    caxs[0].axis('off')

    # VORTICITY
    wmin, wmax = np.amin(w[~mask_hard]), np.amax(w[~mask_hard])
    wmax = max(np.abs(wmin), np.abs(wmax))
    levels = np.linspace(wmin, wmax, 20)

    ct = axs[1].contour(x, y, w, levels=levels, cmap=cmap2)
    cbar = fig.colorbar(ct, cax=caxs[1], ticks=np.round(np.linspace(levels[0], levels[-1], 7)))

    # LABELS
    bbox_dic = dict(boxstyle="round", fc="wheat", ec="none", alpha=0.85)
    axs[0].set_title(r"Streamfunction $\psi$", fontsize=ftSz2)
    axs[1].set_title(r"Vorticity $\omega$", fontsize=ftSz2)
    axs[-1].set_xlabel(r"$x / H_{box}$", fontsize=ftSz2)
    for ax in axs:
        box = Rectangle((din, dbot), lbox, 1., fc="lightgrey", ec="none", alpha=1., zorder=5)
        ax.add_patch(box)
        ax.set_ylabel(r"$y / H_{{box}}$", fontsize=ftSz2)
        ax.axis([0., L, 0., H])

    plt.subplots_adjust(hspace=0.1)
    plt.show()
    return


def plot_max_RE():
    reh = np.empty(nt + 1)
    rew = np.empty(nt + 1)

    for t in tqdm(range(nt+1)):
        p, u, v, w, p_masked, w_masked = read_block(sim)
        u = (u[:, 1:] + u[:, :-1]) / 2.
        v = (v[1:, :] + v[:-1, :]) / 2.
        reh[t] = np.amax(np.abs(u[~mask_corner.T]) + np.abs(v[~mask_corner.T])) * h * RE
        rew[t] = np.amax(np.abs(w[~mask_corner.T])) * h * h * RE

    fig, ax = plt.subplots(1, 1, figsize=(10., 6.))
    fig.tight_layout()
    ax.plot(sim.t, reh, marker=".", markersize=5, color='C0', label='RE_h')
    ax.plot(sim.t, rew, marker=".", markersize=5, color='C1', label='RE_w')

    ax.set_xlabel(r"$t U_{\infty} / H_{box}$", fontsize=ftSz2)
    ax.set_ylabel(r"$Re$", fontsize=ftSz2)
    ax.legend(fontsize=ftSz3)
    ax.grid(ls=':')
    plt.show()

    return


def compute_drag_lift(p, u, v):
    drag_p, drag_f, lift_p, lift_f = 0., 0., 0., 0.

    drag_p += np.trapz(p[n*din-1, n*dbot:n*(dbot+1)], dx=h) - np.trapz(p[n*(din+lbox), n*dbot:n*(dbot+1)], dx=h)
    drag_f += np.trapz(u[n*din:n*(din+lbox)+1, n*dbot] - u[n*din:n*(din+lbox)+1, n*dbot+1])
    drag_f += np.trapz(u[n*din:n*(din+lbox)+1, n*(dbot+1)+1] - u[n*din:n*(din+lbox)+1, n*(dbot+1)])

    lift_p += np.trapz(p[n*din:n*(din+lbox), n*dbot-1], dx=h) - np.trapz(p[n*din:n*(din+lbox), n*(dbot+1)], dx=h)
    lift_f += np.trapz(v[din*n, n*dbot:n*(dbot+1)+1] - v[din*n+1, n*dbot:n*(dbot+1)+1])
    lift_f += np.trapz(v[(din+lbox)*n+1, n*dbot:n*(dbot+1)+1] - v[(din+lbox)*n, n*dbot:n*(dbot+1)+1])

    return drag_p, drag_f, lift_p, lift_f

def plot_drag_lift():
    drag_p, drag_f = np.zeros(nt+1), np.zeros(nt+1)
    lift_p, lift_f = np.zeros(nt+1), np.zeros(nt+1)

    for idx in tqdm(range(nt//3 + 1)):
        p, u, v, _, _, _ = read_block(sim)
        drag1, drag2, lift1, lift2 = compute_drag_lift(p, u, v)
        drag_p[idx], drag_f[idx] = drag1, drag2
        lift_p[idx], lift_f[idx] = lift1, lift2

    fig, axs = plt.subplots(2, 1, figsize=(10., 8.))
    fig.tight_layout()
    axs[0].stackplot(sim.t, drag_p, drag_f, labels=['Pressure', 'Friction'])
    axs[1].stackplot(sim.t, lift_p, lift_f, labels=['Pressure', 'Friction'])

    axs[-1].set_xlabel(r"$t U_{\infty} / H_{box}$", fontsize=ftSz2)
    axs[0].set_ylabel(r"$Drag \;/\; (\rho U_{{\infty}}^2 H_{{box}}^2) $", fontsize=ftSz2)
    axs[1].set_ylabel(r"$Lift \;/\; (\rho U_{{\infty}}^2 H_{{box}}^2) $", fontsize=ftSz2)
    for ax in axs:
        ax.legend(fontsize=ftSz3)
        ax.grid(ls=':')
    
    plt.show()
    return

# kg m / s^2
# F / (rho U0^2 Hbox^2)


if __name__ == "__main__":

    cmap1 = plt.get_cmap("coolwarm")
    cmap2 = plt.get_cmap("bwr")
    cmap3 = plt.get_cmap("Spectral_r")

    path_dir = "../results"
    path_anim = "../anim"

    filename_params = f"{path_dir}/simu_params.txt"
    filename_p = f"{path_dir}/simu_p.txt"
    filename_u = f"{path_dir}/simu_u.txt"
    filename_v = f"{path_dir}/simu_v.txt"
    
    sim = Simulation([filename_params, filename_p, filename_u, filename_v])

    nt, nx, ny, n = sim.nt, sim.nx, sim.ny, sim.n
    T, h, L, H, din, dbot, lbox, RE = sim.T, sim.h, sim.L, sim.H, sim.din, sim.dbot, sim.lbox, sim.RE
    x, y, xx, yy, delta_x, delta_y = sim.x, sim.y, sim.xx, sim.yy, sim.delta_x, sim.delta_y


    dd = 10*h
    # mask_middle = (din < xxm) * (xxm < din+lbox) * (dbot < yym) * (yym < dbot+1)
    mask_corner = (din < xx) * (xx < din + lbox) * (dbot < yy) * (yy < dbot + 1.)
    mask_hard = (din - dd < xx) * (xx < din + lbox + dd) * (dbot - dd < yy) * (yy < dbot + 1. + dd)

    # plot_vorticity([18., 19., 20.], cmap=cmap1)
    # plot_streamlines([18., 19., 20.], cmap=cmap2)
    # plot_average_flow(10., cmap2, cmap1)
    
    # plot_max_RE()
    plot_drag_lift()
