import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

ftSz1, ftSz2, ftSz3 = 20, 17, 15
plt.rcParams['font.family'] = 'monospace'


def read_file(filename):
    with open(filename, 'r') as file:
        scheme = file.readline().strip()
        c, sigma, U_max, L, dt, a = [float(x) for x in file.readline().split()]
    matrix = np.loadtxt(filename, skiprows=2)
    t = matrix[:, 0]
    u = matrix[:, 1:]
    M, N = u.shape
    h = L / N
    return scheme, c, sigma, U_max, L, dt, a, t, u, M, N, h


def plot_soluce(save=False):
    fig, axs = plt.subplots(3, 2, figsize=(9, 10), constrained_layout=True, sharex='all', sharey='all')
    scheme_list = ["E2", "E4", "I4"]
    # fig, axs = plt.subplots(2, 2, figsize=(14, 6.5), constrained_layout=True, sharex='all', sharey='all')
    # scheme_list = ["E6", "I6"]
    
    N_list = [32, 64, 128]
    T_list = [0.500, 1.00]
    prefix = "solution"
    L = 1.

    for i, scheme in enumerate(scheme_list):
        
        for j, (N, mk) in enumerate(zip(N_list, ['.', '.', '.'])):
            scheme, c, sigma, U_max, L, dt, a, t, u, M, N, h = read_file(f"./data/{prefix}_{scheme}_{N}.txt")

            x = np.linspace(-L/2., 3*L/2. - h, 2*N)

            T_idx = [np.argmin(np.abs(t - t_wanted)) for t_wanted in T_list]
            axs[i, 0].plot(x, np.r_[u[T_idx[0]], u[T_idx[0]]], ls='-', marker=mk, color=f'C{j}', label=f'N = {N}', zorder=len(N_list)-j)
            axs[i, 1].plot(x, np.r_[u[T_idx[1]], u[T_idx[1]]], ls='-', marker=mk, color=f'C{j}', zorder=len(N_list)-j)

            # print(f"{scheme:2s} - N = {N:3d}  -->  0.5 =? {T_idx[0] * dt:.3f}  1. =? {T_idx[1] * dt:.3f}")

        n_plot = 500
        x_plot = np.linspace(-L / 2., 3 * L / 2., n_plot)
        f = lambda x_, t_: U_max * np.exp(-np.power((np.fmod((x_ - c * t_ - 3 * L / 2), L) + L / 2) / sigma, 2))

        axs[i, 0].plot(x_plot, f(x_plot, T_list[0]), color='grey', alpha=0.5, lw=5, zorder=0, label='Analytic solution')
        axs[i, 1].plot(x_plot, f(x_plot, T_list[1]), color='grey', alpha=0.5, lw=5, zorder=0, label='Analytic solution')

    lines_labels = [axs[0, 0].get_legend_handles_labels()]
    lines, labels = [sum(lol, []) for lol in zip(*lines_labels)]
    lgd = fig.legend(lines, labels, labelspacing=2.5, bbox_to_anchor=(0.07, -0.99, 0.87, 1.), mode='expand',
                     ncol=4, facecolor='wheat', framealpha=0.25, fancybox=True, fontsize=ftSz3)

    axs[0, 0].set_ylim(-0.15, 1.15)
    for ax in axs.flatten():
        ax.grid(ls=':')
    for T, ax in zip(T_list, axs[0, :]):
        ax.set_title(r"$ct/L \; = \;{:.1f}$".format(T), fontsize=ftSz2)
    for ax in axs[-1, :]:
        ax.set_xlabel("$x/L$", fontsize=ftSz2)
    for scheme, ax in zip(scheme_list, axs[:, 0]):
        ax.set_ylabel(f"$u(x,t)$ - {scheme}", fontsize=ftSz2)
    
    if save:
        fig.savefig("./figures/overview.svg", format="svg", bbox_extra_artists=(lgd,), bbox_inches='tight')
    else:
        plt.show()


def plot_wavepacket(save=False):
    fig, axs = plt.subplots(3, 2, figsize=(10, 10), constrained_layout=True, sharex='all', sharey='all')
    scheme_list = ["E2", "E4", "I4"]
 
    T_list = [0.50, 1.00]
    prefix = "wavepacket"
    N = 128

    for i, scheme in enumerate(scheme_list):
        
        scheme, c, sigma, U_max, L, dt, a, t, u, M, N, h = read_file(f"./data/{prefix}_{scheme}_{N}.txt")
        x = np.linspace(-L/2., 3*L/2. - h, 2*N)

        T_idx = [np.argmin(np.abs(t - t_wanted)) for t_wanted in T_list]
        axs[i, 0].plot(x, np.r_[u[T_idx[0]], u[T_idx[0]]], ls='-', marker='.', label=f'Numerical solution')  # , color=f'C{j}',
        axs[i, 1].plot(x, np.r_[u[T_idx[1]], u[T_idx[1]]], ls='-', marker='.')

        n_plot = 500
        x_plot = np.linspace(-L / 2., 3 * L / 2., n_plot)
        func = lambda x_, t_: U_max * np.exp(-np.power((np.fmod((x_ - c * t_ - 3 * L / 2), L) + L / 2) / sigma, 2))
        f = lambda x_, t_: np.cos(2 * np.pi * 16 * (x_ - c * t_) / L) * func(x_, t_)

        axs[i, 0].plot(x_plot, f(x_plot, T_list[0]), color='grey', alpha=0.5, lw=2, zorder=0, label='Analytic solution')
        axs[i, 1].plot(x_plot, f(x_plot, T_list[1]), color='grey', alpha=0.5, lw=2, zorder=0, label='Analytic solution')

    lines_labels = [axs[0, 0].get_legend_handles_labels()]
    lines, labels = [sum(lol, []) for lol in zip(*lines_labels)]
    lgd = fig.legend(lines, labels, labelspacing=2.5, bbox_to_anchor=(0.225, -0.99, 0.60, 1.), mode='expand',
                     ncol=2, facecolor='wheat', framealpha=0.25, fancybox=True, fontsize=ftSz3)

    axs[0, 0].set_ylim(-0.85, 1.1)
    for ax in axs.flatten():
        ax.grid(ls=':')
    for T, ax in zip(T_list, axs[0, :]):
        ax.set_title(r"$ct/L \; = \;{:.1f}$".format(T), fontsize=ftSz2)
    for ax in axs[-1, :]:
        ax.set_xlabel("$x/L$", fontsize=ftSz2)
    for scheme, ax in zip(scheme_list, axs[:, 0]):
        ax.set_ylabel(f"$u(x,t)$ - {scheme}", fontsize=ftSz2)
    
    if save:
        fig.savefig("./figures/wavepacket.svg", format="svg", bbox_extra_artists=(lgd,), bbox_inches='tight')
    else:
        plt.show()


def plot_problem(save=False):
    fig, axs = plt.subplots(3, 1, figsize=(10, 8), constrained_layout=True, sharex='all', sharey='all')
    T_list = np.array([0.66, 0.682, 0.705])
    scheme = 'E6'
    scheme, c, sigma, U_max, L, dt, a, t, u, M, N, h = read_file(f"./data/solution_{scheme}_128.txt")
    x_loc_packet = np.array([0.150, 0.100, 0.050])

    for i, (ax, t_wanted) in enumerate(zip(axs, T_list)):

        x = np.linspace(-L/2., L/2. - h, N)
        T_idx = np.argmin(np.abs(t - t_wanted))
        ax.plot(x, u[T_idx], ls='-', marker='.', label=f'Numerical solution')  # , color=f'C{j}',

        n_plot = 500
        x_plot = np.linspace(-L / 2., L / 2., n_plot)
        f = lambda x_, t_: U_max * np.exp(-np.power((np.fmod((x_ - c * t_ - 3 * L / 2), L) + L / 2) / sigma, 2))
        ax.plot(x_plot, f(x_plot, t[T_idx]), color='grey', alpha=0.5, lw=2, zorder=0, label='Analytic solution')

        x_loc_regular = c * t[T_idx] + sigma * np.sqrt(np.log(U_max / 0.10)) - L
        ax.vlines(x_loc_regular, -0.10, 0.20, color='C1', lw=3, alpha=0.65, label='Gaussian indicator')
        ax.vlines(x_loc_packet[i], -0.10, 0.20, color='C2', lw=3, alpha=0.65, label='Pulse indicator')

        ax.grid(ls=':')
        ax.set_ylabel(r"$u(x, t={:.3f})$".format(t[T_idx]), fontsize=ftSz2)
        ax.set_xlim([-0.26, 0.19])
        ax.set_ylim([-0.03, 0.16])

    axs[-1].set_xlabel(r"$x\:/\:L$", fontsize=ftSz2)
    axs[0].legend(fontsize=ftSz3, loc='upper center')

    if save:
        fig.savefig("./figures/problem.svg", format="svg", bbox_inches='tight')
    else:
        plt.show()


def plot_soluce_nonuniform(save=False):
    fig, axs = plt.subplots(2, 2, figsize=(10, 8), constrained_layout=True, sharex='all', sharey='row')
    scheme_list = ["E2", "E4", "I4"]  # , "E6", "I6"]
    T_list = [0.50, 1.00]
    L = 1.  # but overwritten

    for i, scheme in enumerate(scheme_list):

        scheme, c, sigma, U_max, L, dt, a, t, v, M, N, h = read_file(f"./data/nonuniform_{scheme}_128.txt")
        xi = np.linspace(-L/2., L/2. - h, N)
        dg = 1 - a * np.cos(2 * np.pi * xi / L)  # x here is xi
        x = xi - a * L / (2 * np.pi) * np.sin(2 * np.pi * xi / L)  # map to physical x
        u = v / dg  # change of variable from v to u

        xi = np.r_[xi, L + xi]
        x = np.r_[x, L + x]
        T_idx = [np.argmin(np.abs(t - t_wanted)) for t_wanted in T_list]
        zorder = len(scheme_list) - i
        axs[0, 0].plot(xi, np.r_[v[T_idx[0]], v[T_idx[0]]], ls='-', marker='', color=f'C{i}', zorder=zorder)
        axs[0, 1].plot(xi, np.r_[v[T_idx[1]], v[T_idx[1]]], ls='-', marker='', color=f'C{i}', zorder=zorder)
        axs[1, 0].plot(x, np.r_[u[T_idx[0]], u[T_idx[0]]], ls='-', marker='', color=f'C{i}', label=f'Scheme {scheme}', zorder=zorder)
        axs[1, 1].plot(x, np.r_[u[T_idx[1]], u[T_idx[1]]], ls='-', marker='', color=f'C{i}', zorder=zorder)

    n_plot = 500
    x_plot = np.linspace(-L / 2., 3 * L / 2., n_plot)
    f = lambda x_, t_: U_max * np.exp(-np.power((np.fmod(np.abs(x_ - c * t_ + L / 2), L) - L / 2) / sigma, 2))
    axs[0, 0].hlines(1.5, -0.5 * L, 1.5 * L, ls='--', color='lightgrey')
    axs[0, 1].hlines(0.5, -0.5 * L, 1.5 * L, ls='--', color='lightgrey')
    axs[1, 0].plot(x_plot, f(x_plot, T_list[0]), color='grey', alpha=0.5, lw=5, zorder=0, label='Analytic solution')
    axs[1, 1].plot(x_plot, f(x_plot, T_list[1]), color='grey', alpha=0.5, lw=5, zorder=0)

    lines_labels = [axs[1, 0].get_legend_handles_labels()]
    lines, labels = [sum(lol, []) for lol in zip(*lines_labels)]
    lgd = fig.legend(lines, labels, labelspacing=2.5, bbox_to_anchor=(0.07, -0.99, 0.87, 1.), mode='expand',
                     ncol=4, facecolor='wheat', framealpha=0.25, fancybox=True, fontsize=ftSz3)

    for ax in axs.flatten():
        ax.grid(ls=':')
    for T, ax in zip(T_list, axs[0, :]):
        ax.set_title(r"$ct/L \; = \;{:.1f}$".format(T), fontsize=ftSz2)
    for ax in axs[0, :]:
        ax.set_xlabel(r"$\xi\:/\:L$", fontsize=ftSz2)
    for ax in axs[1, :]:
        ax.set_xlabel(r"$x\:/\:L$", fontsize=ftSz2)
    axs[0, 0].set_ylabel(r"$v(\xi, t)$", fontsize=ftSz2)
    axs[1, 0].set_ylabel(r"$u(x, t)$", fontsize=ftSz2)

    if save:
        fig.savefig("./figures/nonuniform.svg", format="svg", bbox_extra_artists=(lgd,), bbox_inches='tight')
    else:
        plt.show()


def plot_diagnostic(save=False):
    fig, axs = plt.subplots(3, 3, figsize=(10, 8), constrained_layout=True, sharex='all', sharey='row')
    scheme_list = ["E2", "E4", "I4"]  # , "E6", "I6"]
    
    # fig, axs = plt.subplots(3, 2, figsize=(9, 8), constrained_layout=True, sharex='all', sharey='row')
    # scheme_list = ["E6", "I6"]
    
    N_list = [32, 64, 128]
    for i, scheme in enumerate(scheme_list):
        
        for j, N in enumerate(N_list):
            
            scheme, c, sigma, U_max, L, dt, a, t, u, M, N, h = read_file(f"./data/solution_{scheme}_{N}.txt")
            x = np.linspace(-L/2., L/2. - h, N)
            
            f = lambda x_, t_: U_max * np.exp(-np.power((np.fmod(np.abs(x_ - c * t_ + L / 2), L) - L / 2) / sigma, 2))
            xx, tt = np.meshgrid(x, t)

            I = h / (sigma * U_max) * np.sum(u, axis=1)
            E = h / (sigma * U_max * U_max) * np.sum(u * u, axis=1)
            R = h / (sigma * U_max * U_max) * np.sum(np.power((u - f(xx, tt)), 2), axis=1)

            axs[0, i].plot(t, I, ls='-', marker='.', color=f'C{j}', label=f'N = {N}', zorder=len(N_list)-j)
            axs[1, i].plot(t, E, ls='-', marker='.', color=f'C{j}', zorder=len(N_list)-j)
            axs[2, i].plot(t, R, ls='-', marker='.', color=f'C{j}', zorder=len(N_list)-j)

    lines_labels = [axs[0, 0].get_legend_handles_labels()]
    lines, labels = [sum(lol, []) for lol in zip(*lines_labels)]
    lgd = fig.legend(lines, labels, labelspacing=2.5, bbox_to_anchor=(0.25, -0.99, 0.55, 1.), mode='expand',
                     ncol=3, facecolor='wheat', framealpha=0.25, fancybox=True, fontsize=ftSz3)

    axs[0, 0].set_ylim(1.625, 1.925)
    for ax in axs.flatten():
        ax.grid(ls=':')
    for s, ax in zip(scheme_list, axs[0, :]):
        ax.set_title(f"{s}", fontsize=ftSz2)
    for ax in axs[-1, :]:
        ax.set_xlabel("$ct \: / \: L$", fontsize=ftSz2)
    for scheme, ax in zip([r"$I_h^n$", r"$E_h^n$", r"$R_h^n$"], axs[:, 0]):
        ax.set_ylabel(f"{scheme}", fontsize=ftSz2)
    
    if save:
        fig.savefig("./figures/diagnostic.svg", format="svg", bbox_extra_artists=(lgd,), bbox_inches='tight')
    else:
        plt.show()


def order_convergence(save=False):
    fig, ax = plt.subplots(1, 1, figsize=(9., 5.), constrained_layout=True)
    scheme_list = ["E2", "E4", "I4"]  # , "E6", "I6"]
    alpha_list = [1., 1., 1., 0.5, 0.5]
    N_list = [32, 64, 128]
    t_wanted = 0.5
    h_list = []  # but overwritten

    for i, scheme in enumerate(scheme_list):
        
        R_list = np.zeros(3)
        h_list = np.zeros(3)
        for j, N in enumerate(N_list):
            
            scheme, c, sigma, U_max, L, dt, a, t, u, M, N, h = read_file(f"./data/solution_{scheme}_{N}.txt")
            x = np.linspace(-L/2., L/2. - h, N)
            
            f = lambda x_, t_: U_max * np.exp(-np.power((np.fmod(np.abs(x_ - c * t_ + L / 2), L) - L / 2) / sigma, 2))
            
            T_idx = np.argmin(np.abs(t - t_wanted))
            h_list[j] = h / sigma
            R_list[j] = h / (sigma * U_max * U_max) * np.sum(np.power((u[T_idx] - f(x, t[T_idx])), 2))

        ax.loglog(h_list, R_list, ls='-', marker='o', color=f'C{i}', label=f'{scheme}', alpha=alpha_list[i])
    
    ax.loglog(h_list, h_list ** 2, ls='--', color='black', label='order $h^2$')
    ax.loglog(h_list, h_list ** 4, ls='-.', color='black', label='order $h^4$')

    ax.grid(ls=':', which="both")
    ax.set_xlabel(r"$h\:/\:\sigma$", fontsize=ftSz2)
    ax.set_ylabel(r"Global error at $ct\:/\:L \:=\: {:.1f}$".format(t_wanted), fontsize=ftSz2)
    ax.legend(fontsize=ftSz3)
    
    if save:
        fig.savefig("./figures/order_small_fig.svg", format="svg", bbox_inches='tight')
    else:
        plt.show()


def animation_soluce(blit=False):
    def init():
        exact.set_data(x_plot, f(x_plot, 0))
        time_text.set_text(time_template.format(0))
        line.set_data(x, u[0, :])

        return tuple([line, exact, time_text])

    def animate(t_idx):
        exact.set_data(x_plot, f(x_plot, t[t_idx]))
        time_text.set_text(time_template.format(t[t_idx]))
        line.set_data(x, u[t_idx, :])

        return tuple([line, exact, time_text])

    scheme, c, sigma, U_max, L, dt, a, t, u, M, N, h = read_file("./data/solution.txt")
    x = np.linspace(-L/2., L/2. - h, N)
    dg = 1 - a * np.cos(2 * np.pi * x / L)  # x here is xi
    x = x - a * L / (2 * np.pi) * np.sin(2 * np.pi * x / L)  # map to physical x
    u /= dg  # change of variable from v to u

    func = lambda x_, t_: U_max * np.exp(-np.power((np.fmod((x_ - c * t_ - L / 2), L) + L / 2) / sigma, 2))
    if np.any(u[0] < 0):
        f = lambda x_, t_: np.cos(2 * np.pi * 16 * (x_ - c * t_) / L) * func(x_, t_)
    else:
        f = func

    n_plot = 200
    x_plot = np.linspace(-L / 2., L / 2., n_plot)

    fig, ax = plt.subplots(1, 1, figsize=(10, 6), constrained_layout=True, num='Animation of the solution')

    time_template = r'$t = {:.2f} \;[s]$'
    time_text = ax.text(0.03, 0.88, '', fontsize=17, transform=ax.transAxes)
    ax.text(0.03, 0.94, 'FD scheme: {:s}'.format(scheme), fontsize=17, transform=ax.transAxes)

    line, = ax.plot(x, u[0, :], ls='-', marker='.', color='C0', label='Numerical solution')
    exact, = ax.plot(x_plot, f(x_plot, 0), color='C1', alpha=0.5, lw=5, zorder=0, label='Analytic solution')
    ax.legend(fontsize=ftSz2, loc='upper right')

    ax.grid(ls=':')
    ax.set_xlabel("x", fontsize=ftSz2)
    ax.set_ylabel("u(x,t)", fontsize=ftSz2)

    # to animate
    _ = FuncAnimation(fig, animate, M, interval=20, blit=blit, init_func=init, repeat_delay=3000)

    # to get only one frame at t = i
    # i = M-1 ; init() ; animate(i)
    plt.show()


if __name__ == "__main__":
    
    save_global = False
    plt.rcParams["text.usetex"] = save_global

    print("{:15s}".format("Anim running"), end="\r")
    animation_soluce(blit=False)
    print("{:15s}".format("Plot 1 / 6"), end="\r")
    plot_soluce(save_global)
    print("{:15s}".format("Plot 2 / 6"), end="\r")
    plot_diagnostic(save_global)
    print("{:15s}".format("Plot 3 / 6"), end="\r")
    order_convergence(save_global)
    print("{:15s}".format("Plot 4 / 6"), end="\r")
    plot_soluce_nonuniform(save_global)
    print("{:15s}".format("Plot 5 / 6"), end="\r")
    plot_wavepacket(save_global)
    print("{:15s}".format("Plot 6 / 6"), end="\r")
    plot_problem(save_global)
    print("{:15s}".format("Job done"))

