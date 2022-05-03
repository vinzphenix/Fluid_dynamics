# import os
# sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), 'run.py')))

from unittest import case
from run_new import *


def plot_vorticity(t_list, cmap, save=False):

    n_plots = len(t_list)
    indices = [np.argmin(np.abs(this_t - sim.t)) for this_t in t_list]
        
    fig, axs = plt.subplots(n_plots, 1, figsize=(10., 10.), sharex="all", sharey="all")
    fig.tight_layout()
    caxs = [make_colorbar_with_padding(axs[i]) for i in range(n_plots)]
    
    i = 0
    for t in tqdm(range(nt+1)):
        w, = read_block(sim, needed={"w": True})
        if t not in indices:
            continue

        w = np.ma.masked_array(w.T, mask_corner)
        
        wmin, wmax = np.amin(w[~mask_hard]), np.amax(w[~mask_hard])
        wmax = np.round(max(np.abs(wmin), np.abs(wmax)))
        wmin = -wmax
        wmin, wmax = -25, 25

        scaling, nStreamLines = 0.5, 40
        wavg, wdif = (wmax + wmin) / 2., wmax - wmin
        levels = np.linspace(wmin + 0.01 * wdif, wmax - 0.01 * wdif, nStreamLines)
        levels = levels - scaling * wdif / (2. * np.pi) * np.sin(2. * np.pi / wdif * (levels - wavg))
        
        # levels = np.linspace(wmax * 0.05, np.round(wmax * 0.95, 0), 10)
        # levels = np.r_[-levels[::-1], levels]
        # levels = np.linspace(-wmax, wmax, 20)
        # levels = levels - 0.5 * (2*wmax) / (2. * np.pi) * np.sin(2. * np.pi / (2 * wmax) * (levels))

        ct = axs[i].contour(x + delta_x[t], y + delta_y[t], w, linewidths=2., levels=levels, cmap=cmap)
        cbar = fig.colorbar(ct, cax=caxs[i], spacing = 'proportional', ticks=np.round(np.linspace(levels[0], levels[-1], 7)))
        i += 1
        if i == n_plots:
            break

    bbox_dic = dict(boxstyle="round", fc="wheat", ec="none", alpha=0.85)
    axs[-1].set_xlabel(r"$x / H_{box}$", fontsize=ftSz2)
    for idx, ax in zip(indices, axs):
        box = Rectangle((din + delta_x[idx], dbot + delta_y[idx]), lbox, 1., fc="grey", ec="none", alpha=1., zorder=5)
        print(din + delta_x[idx])
        ax.add_patch(box)
        ax.text(0.03,0.85, r"$t = {:5.2f}$".format(sim.t[idx]), fontsize=ftSz2, transform=ax.transAxes, bbox=bbox_dic)
        ax.set_ylabel(r"$y / H_{{box}}$", fontsize=ftSz2)
        ax.axis([0., L, 0., H])
        ax.set_aspect("equal")

    if save:
        fig.savefig(f"../figures/vorticity_{case_name}.svg", format="svg", bbox_inches='tight')
    else:
        plt.show()


def plot_streamlines(t_list, cmap, save=False):

    n_plots = len(t_list)
    indices = [np.argmin(np.abs(this_t - sim.t)) for this_t in t_list]
    x, y, delta_x, delta_y = sim.x, sim.y, sim.delta_x, sim.delta_y

    fig, axs = plt.subplots(n_plots, 1, figsize=(10., 10.), sharex="all", sharey="all")
    fig.tight_layout()
    caxs = [make_colorbar_with_padding(axs[i]) for i in range(n_plots)]

    i = 0
    cbars = []
    for idx in tqdm(range(nt+1)):
        u, v = read_block(sim, needed={"u": True, "v": True})
        if idx not in indices:
            continue
        u = (u[:, 1:] + u[:, :-1]) / 2.
        v = (v[1:, :] + v[:-1, :]) / 2.
        speed = np.hypot(u.T, v.T)
        lw = 2 * np.abs(1. - speed)
        strm = axs[i].streamplot(x + delta_x[idx], y + delta_y[idx], u.T, v.T, density=2., color=speed, cmap=cmap, linewidth=lw)
        cbars.append(fig.colorbar(strm.lines, cax=caxs[i]))
        cbars[i].ax.set_ylabel(r"$\|\: (u,\, v) \:\|$", fontsize=ftSz2)
        i += 1
        if i == n_plots:
            break

    bbox_dic = dict(boxstyle="round", fc="lightgrey", ec="none", alpha=0.85)
    axs[-1].set_xlabel(r"$x / H_{box}$", fontsize=ftSz2)
    for idx, ax in zip(indices, axs):
        box = Rectangle((din + delta_x[idx], dbot + delta_y[idx]), lbox, 1., fc="lightgrey", ec="none", alpha=1., zorder=5)
        ax.add_patch(box)
        ax.text(0.85,0.85, r"$t = {:5.2f}$".format(sim.t[idx]), fontsize=ftSz2, transform=ax.transAxes, bbox=bbox_dic)
        ax.set_ylabel(r"$y / H_{{box}}$", fontsize=ftSz2)
        ax.axis([0., L, 0., H])
        ax.set_aspect("equal")

    if save:
        fig.savefig(f"../figures/streamlines_{case_name}.svg", format="svg", bbox_inches='tight')
    else:
        plt.show()


def plot_average_flow(cmap1, cmap2, compute=False, t_start=20., save=False):
    if sim.T < t_start:
        return

    fig, axs = plt.subplots(2, 1, figsize=(10., 8.), sharex="all", sharey="all")
    fig.tight_layout()
    caxs = [make_colorbar_with_padding(axs[i]) for i in range(2)]

    if compute:
        idx_start = np.argmin(np.abs(sim.t - t_start))

        p_avg = np.zeros((nx, ny))
        u_avg = np.zeros((nx+1, ny+2))
        v_avg = np.zeros((nx+2, ny+1))

        for idx in tqdm(range(nt + 1)):
            p, u, v, w = read_block(sim, needed={"p": True, "u": True, "v": True, "w": True})
            # print(np.amax(np.abs((u[:, 1:] + u[:, :-1]) / 2.) + np.abs((v[1:, :] + v[:-1, :]) / 2.)))
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

    else:  # Using values computed during simulation
        u_avg, v_avg = sim.u_avg, sim.v_avg
        w = ((v_avg[1:, :] - v_avg[:-1, :]) - (u_avg[:, 1:] - u_avg[:, :-1])) / sim.h
        w = np.ma.masked_array(w.T, mask_corner)

    u_avg = (u_avg[:, 1:] + u_avg[:, :-1]) / 2.
    v_avg = (v_avg[1:, :] + v_avg[:-1, :]) / 2.

    # STREAMLINES
    # psi = axs[0].contour(x, y, u_sum.T, 10, cmap=cmap)
    # cbar1 = fig.colorbar(psi, cax=caxs[0])
    speed = np.hypot(u_avg.T, v_avg.T)
    lw = 2 * np.abs(1. - speed)
    strm = axs[0].streamplot(x, y, u_avg.T, v_avg.T, density=2., linewidth=lw, color=speed, cmap=cmap1)
    cbar1 = fig.colorbar(strm.lines, cax=caxs[0])
    cbar1.ax.set_ylabel(r"$\|\: (u,\, v) \:\|$", fontsize=ftSz2)

    # VORTICITY
    wmin, wmax = np.amin(w[~mask_hard]), np.amax(w[~mask_hard])
    wmax = max(np.abs(wmin), np.abs(wmax))
    levels = np.linspace(-0.75 * wmax, 0.75 * wmax, 50)
    ct = axs[1].contour(x, y, w, levels=levels, cmap=cmap2)
    cbar2 = fig.colorbar(ct, cax=caxs[1], ticks=np.round(np.linspace(levels[0], levels[-1], 7)))
    cbar2.ax.set_ylabel(r"$\omega$", fontsize=ftSz2)

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

    if save:
        fig.savefig(f"../figures/avg_flow_{case_name}.svg", format="svg", bbox_inches='tight')
    else:
        plt.show()


def plot_max_RE(compute=False, save=False):
    fig, ax = plt.subplots(1, 1, figsize=(10., 5.))
    fig.tight_layout()

    if compute:
        reh = np.empty(nt + 1)
        rew = np.empty(nt + 1)
        for t in tqdm(range(nt+1)):
            u, v, w = read_block(sim, needed={"u": True, "v": True, "w": True})
            u = (u[:, 1:] + u[:, :-1]) / 2.
            v = (v[1:, :] + v[:-1, :]) / 2.
            reh[t] = np.amax(np.abs(u[~mask_corner.T]) + np.abs(v[~mask_corner.T])) * h * RE
            rew[t] = np.amax(np.abs(w[~mask_corner.T])) * h * h * RE

        ax.plot(sim.t, reh, marker=".", markersize=5, color='C0', label='$RE_h$')
        ax.plot(sim.t, rew, marker=".", markersize=5, color='C1', label='$RE_{\omega}$')

    t_start = 0.1
    idx_st = np.argmin(np.abs(sim.t_simu - t_start))
    t_slice = sim.t_simu[idx_st:]
    line1, = ax.plot(t_slice, sim.reh[idx_st:], color='C0', label="$Re_h$")
    line2, = ax.plot(t_slice, sim.rew[idx_st:], color='C1', label="$Re_{\omega}$")

    ax.set_xlabel(r"$t U_{\infty} / H_{box}$", fontsize=ftSz2)
    ax.set_ylabel(r"$Re$", fontsize=ftSz2)
    ymin, ymax = ax.get_ylim()
    line1.set_data(sim.t_simu[1:], sim.reh[1:])
    line2.set_data(sim.t_simu[1:], sim.rew[1:])

    ax.set_ylim([ymin, ymax])
    ax.legend(fontsize=ftSz3)
    ax.grid(ls=':')

    if save:
        fig.savefig(f"../figures/mesh_Re_{case_name}.svg", format="svg", bbox_inches='tight')
    else:
        plt.show()


def compute_drag_lift(p, u, v):
    drag_p, drag_f, lift_p, lift_f = 0., 0., 0., 0.

    drag_p += np.trapz(p[n*din-1, n*dbot:n*(dbot+1)], dx=h) - np.trapz(p[n*(din+lbox), n*dbot:n*(dbot+1)], dx=h)
    drag_f += np.trapz(u[n*din:n*(din+lbox)+1, n*dbot] - u[n*din:n*(din+lbox)+1, n*dbot+1]) / sim.RE
    drag_f += np.trapz(u[n*din:n*(din+lbox)+1, n*(dbot+1)+1] - u[n*din:n*(din+lbox)+1, n*(dbot+1)]) / sim.RE

    lift_p += np.trapz(p[n*din:n*(din+lbox), n*dbot-1], dx=h) - np.trapz(p[n*din:n*(din+lbox), n*(dbot+1)], dx=h)
    lift_f += np.trapz(v[din*n, n*dbot:n*(dbot+1)+1] - v[din*n+1, n*dbot:n*(dbot+1)+1]) / sim.RE
    lift_f += np.trapz(v[(din+lbox)*n+1, n*dbot:n*(dbot+1)+1] - v[(din+lbox)*n, n*dbot:n*(dbot+1)+1]) / sim.RE

    return drag_p, drag_f, lift_p, lift_f

def plot_drag_lift(compute=False, save=False):
    fig, axs = plt.subplots(2, 1, figsize=(10., 6.))
    fig.tight_layout()

    # Using saved fields (not as good, because times sampled)
    if compute:
        drag_p, drag_f = np.zeros(nt+1), np.zeros(nt+1)
        lift_p, lift_f = np.zeros(nt+1), np.zeros(nt+1)
        
        for idx in tqdm(range(nt + 1)):
            p, u, v = read_block(sim, needed={"p": True, "u": True, "v": True})
            drag1, drag2, lift1, lift2 = compute_drag_lift(p, u, v)
            drag_p[idx], drag_f[idx] = 2*drag1, 2*drag2
            lift_p[idx], lift_f[idx] = 2*lift1, 2*lift2
        
        axs[0].stackplot(sim.t, drag_p, drag_f, labels=['Pressure', 'Friction'])
        axs[1].stackplot(sim.t, lift_p, lift_f, labels=['Pressure', 'Friction'])
        axs[0].plot(sim.t, drag_p, 'o', color='C0', label='Pressure')
        axs[0].plot(sim.t, drag_f, 'o', color='C1', label='Friction')
        axs[0].plot(sim.t, drag_p+drag_f, 'o', color='black', label='Total')
        axs[1].plot(sim.t, lift_p, 'o', color='C0', label='Pressure')
        axs[1].plot(sim.t, lift_f, 'o', color='C1', label='Friction')
        axs[1].plot(sim.t, lift_p+lift_f, 'o', color='black', label='Total')

    t_ref1, t_ref2 = 0.2, 50.
    st = np.argmin(np.abs(sim.t_simu - t_ref1))
    fn = np.argmin(np.abs(sim.t_simu - t_ref2))
    t_full = sim.t_simu[st:fn]

    drag_tot = sim.drag_p + sim.drag_f
    axs[0].plot(t_full, sim.drag_p[st:fn], color='C2')  #, label='Pressure')
    axs[0].plot(t_full, sim.drag_f[st:fn], color='C1')  #, label='Friction')
    axs[0].plot(t_full, drag_tot[st:fn], color='C0')  #, label='Total')
    lift_tot = sim.lift_p + sim.lift_f
    axs[1].plot(t_full, sim.lift_p[st:fn], color='C2', label='Pressure')
    axs[1].plot(t_full, sim.lift_f[st:fn], color='C1', label='Friction')
    axs[1].plot(t_full, lift_tot[st:fn], color='C0', label='Total')

    # axs[0].stackplot(t_full, sim.drag_p[st:], sim.drag_f[st:], labels=['Pressure', 'Friction'])
    # axs[1].stackplot(t_full, sim.lift_p[st:], sim.lift_f[st:], labels=['Pressure', 'Friction'])

    axs[-1].set_xlabel(r"$t U_{\infty} / H_{box}$", fontsize=ftSz2)
    axs[0].set_ylabel(r"$C_d$", fontsize=ftSz2)
    axs[1].set_ylabel(r"$C_l$", fontsize=ftSz2)
    # axs[0].set_ylabel(r"$Drag \;/\; (\rho U_{{\infty}}^2 H_{{box}}^2) $", fontsize=ftSz2)
    # axs[1].set_ylabel(r"$Lift \;/\; (\rho U_{{\infty}}^2 H_{{box}}^2) $", fontsize=ftSz2)
    
    axs[1].legend(fontsize=ftSz3, bbox_to_anchor=(0.1, -0.30, 0.8, 0.1), mode="expand", fancybox=True, shadow=True, ncol=3)
    for ax in axs:
        ax.grid(ls=':')
    
    if save:
        fig.savefig(f"../figures/drag_lift_{case_name}.svg", format="svg", bbox_inches='tight')
    else:
        plt.show()

# kg m / s^2
# F / (rho U0^2 Hbox^2)



if __name__ == "__main__":
    save_global = True
    ftSz1, ftSz2, ftSz3 = 25, 20, 17
    plt.rcParams["text.usetex"] = save_global
    plt.rcParams['font.family'] = 'serif'  # monospace

    cmap1 = plt.get_cmap("coolwarm")
    cmap2 = plt.get_cmap("bwr")
    cmap3 = plt.get_cmap("Spectral_r")

    case_list = ["case_1", "case_2", "case_3", "case_4", "hot_box", "hot_cold", "no_slip",
                 "Re_2k_case_1", "Re_2k_eckert", "Re_2k_hot_box", "Re_100_hot_cold", "vertical"]
    
    # for case_name in case_list[:4]:
    for case_name in ["vertical"]:
        path_res = "/media/vincelinux/LaCie/LMECA2660/results/" + case_name
        # path_res = "../results/" + case_name

        filename_params = f"{path_res}/simu_params.txt"
        filename_stats = f"{path_res}/simu_stats.txt"
        filename_p = f"{path_res}/simu_p.txt"
        filename_u = f"{path_res}/simu_u.txt"
        filename_v = f"{path_res}/simu_v.txt"
        filename_T = f"{path_res}/simu_T.txt"
        filename_u_avg = f"{path_res}/simu_u_avg.txt"
        filename_v_avg = f"{path_res}/simu_v_avg.txt"

        sim = Simulation([filename_params, filename_p, filename_u, filename_v, filename_T,
                        filename_stats, filename_u_avg, filename_v_avg], analyze=True)

        nt, nx, ny, n = sim.nt, sim.nx, sim.ny, sim.n
        T, h, L, H, din, dbot, lbox, RE = sim.T, sim.h, sim.L, sim.H, sim.din, sim.dbot, sim.lbox, sim.RE
        x, y, xx, yy, delta_x, delta_y = sim.x, sim.y, sim.xx, sim.yy, sim.delta_x, sim.delta_y


        dd = 3*h
        # mask_middle = (din < xxm) * (xxm < din+lbox) * (dbot < yym) * (yym < dbot+1)
        mask_corner = (din < xx) * (xx < din + lbox) * (dbot < yy) * (yy < dbot + 1.)
        mask_hard = (din - dd < xx) * (xx < din + lbox + dd) * (dbot - dd < yy) * (yy < dbot + 1. + dd)


        # CAN ONLY DO ONE AT A TIME

        # plot_vorticity([12.5e-0, 25.e-0, 50.e-0], cmap=cmap1, save=save_global)
        plot_vorticity([20.e-0, 22.e-0, 24.e-0], cmap=cmap1, save=save_global)  # case2


        # plot_streamlines([10., 20., 30.], cmap=cmap2, save=save_global)  # 0.2, 0.5, 1.
        # plot_average_flow(cmap2, cmap1, compute=False, t_start=20., save=save_global)
        # plot_max_RE(compute=False, save=save_global)
        # plot_drag_lift(compute=False, save=save_global)
