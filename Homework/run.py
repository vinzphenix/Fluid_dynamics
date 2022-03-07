import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

ftSz1, ftSz2, ftSz3 = 20, 15, 12
plt.rcParams['font.family'] = 'monospace'

filename = "./data/solution.txt"


def plot_soluce():
    fig, axs = plt.subplots(3, 2, figsize=(14, 10), constrained_layout=True, sharex='all', sharey='all')
    scheme_list = ["E2", "E4", "I4"]
    # fig, axs = plt.subplots(2, 2, figsize=(14, 6.5), constrained_layout=True, sharex='all', sharey='all')
    # scheme_list = ["E6", "I6"]
    
    N_list = [32, 64, 128]
    T_list = [0.525, 1.006]
    prefix = "solution"

    for i, scheme in enumerate(scheme_list):
        
        for j, N in enumerate(N_list):
            
            with open(f"./data/{prefix}_{scheme}_{N}.txt", 'r') as file:
                scheme = file.readline().strip()
                c, sigma, U_max, L, dt, a = [float(x) for x in file.readline().split()]
                matrix = np.loadtxt(file)

            t = matrix[:, 0]
            u = matrix[:, 1:]
            M, N = u.shape
            h = L / N
            x = np.linspace(-L/2., 3*L/2. - h, 2*N)

            T_idx = [np.argmin(np.abs(t - t_wanted)) for t_wanted in T_list]
            axs[i, 0].plot(x, np.r_[u[T_idx[0]], u[T_idx[0]]], ls='-', marker='.', color=f'C{j}', label=f'N = {N}', zorder=len(N_list)-j)
            axs[i, 1].plot(x, np.r_[u[T_idx[1]], u[T_idx[1]]], ls='-', marker='.', color=f'C{j}', zorder=len(N_list)-j)

            print(f"{scheme:2s} - N = {N:3d}  -->  0.5 =? {T_idx[0] * dt:.3f}  1. =? {T_idx[1] * dt:.3f}")

        n_plot = 500
        x_plot = np.linspace(-L / 2., 3 * L / 2., n_plot)
        f = lambda x_, t_: U_max * np.exp(-np.power((np.fmod((x_ - c * t_ - 3 * L / 2), L) + L / 2) / sigma, 2))

        axs[i, 0].plot(x_plot, f(x_plot, T_list[0]), color='grey', alpha=0.5, lw=5, zorder=0, label='Analytic solution')
        axs[i, 1].plot(x_plot, f(x_plot, T_list[1]), color='grey', alpha=0.5, lw=5, zorder=0, label='Analytic solution')

    lines_labels = [axs[0, 0].get_legend_handles_labels()]
    lines, labels = [sum(lol, []) for lol in zip(*lines_labels)]
    lgd = fig.legend(lines, labels, labelspacing=2.5, bbox_to_anchor=(0.2, -0.99, 0.6, 1.), mode='expand',
                     ncol=4, facecolor='wheat', framealpha=0.25, fancybox=True, fontsize=ftSz2)

    axs[0, 0].set_ylim(-0.1, 1.1)
    for ax in axs.flatten():
        ax.grid(ls=':')
    for T, ax in zip(T_list, axs[0, :]):
        ax.set_title(r"$ct/L \; = \;{:.3f}$".format(T), fontsize=ftSz2)
    for ax in axs[-1, :]:
        ax.set_xlabel("$x/L$", fontsize=ftSz2)
    for scheme, ax in zip(scheme_list, axs[:, 0]):
        ax.set_ylabel(f"u(x,t) - {scheme}", fontsize=ftSz2)
    
    # fig.savefig("./figures/overview_1.svg", format="svg", bbox_extra_artists=(lgd,), bbox_inches='tight')
    plt.show()


def plot_wavepacket():
    fig, axs = plt.subplots(3, 2, figsize=(14, 10), constrained_layout=True, sharex='all', sharey='all')
    scheme_list = ["E2", "E4", "I4"]
 
    T_list = [0.503, 1.006]
    prefix = "wavepacket"
    N = 128

    for i, scheme in enumerate(scheme_list):
        
        with open(f"./data/{prefix}_{scheme}_{N}.txt", 'r') as file:
            scheme = file.readline().strip()
            c, sigma, U_max, L, dt, a = [float(x) for x in file.readline().split()]
            matrix = np.loadtxt(file)

        t = matrix[:, 0]
        u = matrix[:, 1:]
        M, N = u.shape
        h = L / N
        x = np.linspace(-L/2., 3*L/2. - h, 2*N)

        T_idx = [np.argmin(np.abs(t - t_wanted)) for t_wanted in T_list]
        axs[i, 0].plot(x, np.r_[u[T_idx[0]], u[T_idx[0]]], ls='-', marker='.', label=f'Numercial solution')  # , color=f'C{j}',
        axs[i, 1].plot(x, np.r_[u[T_idx[1]], u[T_idx[1]]], ls='-', marker='.')

        n_plot = 500
        x_plot = np.linspace(-L / 2., 3 * L / 2., n_plot)
        func = lambda x_, t_: U_max * np.exp(-np.power((np.fmod((x_ - c * t_ - 3 * L / 2), L) + L / 2) / sigma, 2))
        f = lambda x_, t_: np.cos(2 * np.pi * 16 * (x_ - c * t_) / L) * func(x_, t_) 

        axs[i, 0].plot(x_plot, f(x_plot, T_list[0]), color='grey', alpha=0.5, lw=2, zorder=0, label='Analytic solution')
        axs[i, 1].plot(x_plot, f(x_plot, T_list[1]), color='grey', alpha=0.5, lw=2, zorder=0, label='Analytic solution')

    lines_labels = [axs[0, 0].get_legend_handles_labels()]
    lines, labels = [sum(lol, []) for lol in zip(*lines_labels)]
    lgd = fig.legend(lines, labels, labelspacing=2.5, bbox_to_anchor=(0.2, -0.99, 0.6, 1.), mode='expand',
                     ncol=2, facecolor='wheat', framealpha=0.25, fancybox=True, fontsize=ftSz2)

    axs[0, 0].set_ylim(-0.85, 1.1)
    for ax in axs.flatten():
        ax.grid(ls=':')
    for T, ax in zip(T_list, axs[0, :]):
        ax.set_title(r"$ct/L \; = \;{:.3f}$".format(T), fontsize=ftSz2)
    for ax in axs[-1, :]:
        ax.set_xlabel("$x/L$", fontsize=ftSz2)
    for scheme, ax in zip(scheme_list, axs[:, 0]):
        ax.set_ylabel(f"u(x,t) - {scheme}", fontsize=ftSz2)
    
    # fig.savefig("./figures/wavepacket_1.svg", format="svg", bbox_extra_artists=(lgd,), bbox_inches='tight')
    plt.show()


def plot_soluce_nonuniform():
    fig, axs = plt.subplots(2, 2, figsize=(14, 10), constrained_layout=True, sharex='all', sharey='row')
    scheme_list = ["E2", "E4", "I4"]  # , "E6", "I6"]

    T_list = [0.50312, 1.00078]
    for i, scheme in enumerate(scheme_list):
            
        with open(f"./data/nonuniform_{scheme}_128.txt", 'r') as file:
            scheme = file.readline().strip()
            c, sigma, U_max, L, dt, a = [float(x) for x in file.readline().split()]
            matrix = np.loadtxt(file)

        t = matrix[:, 0]
        v = matrix[:, 1:]
        M, N = v.shape
        h = L / N
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
        axs[1, 0].plot(x, np.r_[u[T_idx[0]], u[T_idx[0]]], ls='-', marker='', color=f'C{i}', label=f'Scheme = {scheme}', zorder=zorder)
        axs[1, 1].plot(x, np.r_[u[T_idx[1]], u[T_idx[1]]], ls='-', marker='', color=f'C{i}', zorder=zorder)

    n_plot = 500
    x_plot = np.linspace(-L / 2., 3 * L / 2., n_plot)
    f = lambda x_, t_: U_max * np.exp(-np.power((np.fmod(np.abs(x_ - c * t_ + L / 2), L) - L / 2) / sigma, 2))
    axs[1, 0].plot(x_plot, f(x_plot, T_list[0]), color='grey', alpha=0.5, lw=5, zorder=0, label='Analytic solution')
    axs[1, 1].plot(x_plot, f(x_plot, T_list[1]), color='grey', alpha=0.5, lw=5, zorder=0)

    lines_labels = [axs[1, 0].get_legend_handles_labels()]
    lines, labels = [sum(lol, []) for lol in zip(*lines_labels)]
    lgd = fig.legend(lines, labels, labelspacing=2.5, bbox_to_anchor=(0.1, -0.99, 0.8, 1.), mode='expand',
                     ncol=4, facecolor='wheat', framealpha=0.25, fancybox=True, fontsize=ftSz2)

    for ax in axs.flatten():
        ax.grid(ls=':')
    for T, ax in zip(T_list, axs[0, :]):
        ax.set_title(r"$ct/L \; = \;{:.3f}$".format(T), fontsize=ftSz2)
    for ax in axs[0, :]:
        ax.set_xlabel(r"$\xi/L$", fontsize=ftSz2)
    for ax in axs[1, :]:
        ax.set_xlabel(r"$x/L$", fontsize=ftSz2)
    axs[0, 0].set_ylabel(r"$v(\xi, t)$", fontsize=ftSz2)
    axs[1, 0].set_ylabel(r"$u(x, t)$", fontsize=ftSz2)

    # fig.savefig("./figures/nonuniform_1.svg", format="svg", bbox_extra_artists=(lgd,), bbox_inches='tight')
    plt.show()


def plot_diagnostic():
    fig, axs = plt.subplots(3, 3, figsize=(12, 8), constrained_layout=True, sharex='all', sharey='row')
    scheme_list = ["E2", "E4", "I4"]  # , "E6", "I6"]
    
    # fig, axs = plt.subplots(3, 2, figsize=(9, 8), constrained_layout=True, sharex='all', sharey='row')
    # scheme_list = ["E6", "I6"]
    
    N_list = [32, 64, 128]
    for i, scheme in enumerate(scheme_list):
        
        for j, N in enumerate(N_list):
            
            with open(f"./data/solution_{scheme}_{N}.txt", 'r') as file:
                scheme = file.readline().strip()
                c, sigma, U_max, L, dt, a = [float(x) for x in file.readline().split()]
                matrix = np.loadtxt(file)

            t = matrix[:, 0]
            u = matrix[:, 1:]
            M, N = u.shape
            h = L / N
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
    lgd = fig.legend(lines, labels, labelspacing=2.5, bbox_to_anchor=(0.2, -0.99, 0.6, 1.), mode='expand',
                     ncol=3, facecolor='wheat', framealpha=0.25, fancybox=True, fontsize=ftSz2)

    axs[0, 0].set_ylim(1.625, 1.925)
    for ax in axs.flatten():
        ax.grid(ls=':')
    for s, ax in zip(scheme_list, axs[0, :]):
        ax.set_title(f"{s}", fontsize=ftSz2)
    for ax in axs[-1, :]:
        ax.set_xlabel("ct/L", fontsize=ftSz2)
    for scheme, ax in zip(["I_h", "E_h", "R_h"], axs[:, 0]):
        ax.set_ylabel(f"{scheme}", fontsize=ftSz2)
    
    # fig.savefig("./figures/diagnostic.svg", format="svg", bbox_extra_artists=(lgd,), bbox_inches='tight')
    plt.show()


def order_convergence():
    fig, ax = plt.subplots(1, 1, figsize=(10., 6.), constrained_layout=True)
    scheme_list = ["E2", "E4", "I4", "E6", "I6"]
    alpha_list = [1., 1., 1., 0.5, 0.5]
    N_list = [32, 64, 128]
    t_wanted = 0.525

    for i, scheme in enumerate(scheme_list):
        
        R_list = np.zeros(3)
        h_list = np.zeros(3)
        for j, N in enumerate(N_list):
            
            with open(f"./data/solution_{scheme}_{N}.txt", 'r') as file:
                scheme = file.readline().strip()
                c, sigma, U_max, L, dt, a = [float(x) for x in file.readline().split()]
                matrix = np.loadtxt(file)

            t = matrix[:, 0]
            u = matrix[:, 1:]
            M, N = u.shape
            h = L / N
            x = np.linspace(-L/2., L/2. - h, N)
            
            f = lambda x_, t_: U_max * np.exp(-np.power((np.fmod(np.abs(x_ - c * t_ + L / 2), L) - L / 2) / sigma, 2))
            
            T_idx = np.argmin(np.abs(t - t_wanted))
            h_list[j] = h / sigma
            R_list[j] = h / (sigma * U_max * U_max) * np.sum(np.power((u[T_idx] - f(x, t[T_idx])), 2))

        ax.loglog(h_list, R_list, ls='-', marker='o', color=f'C{i}', label=f'{scheme}', alpha=alpha_list[i])
    
    ax.loglog(h_list, h_list ** 2, ls='--', color='black', label='order $h^2$')
    ax.loglog(h_list, h_list ** 4, ls='-.', color='black', label='order $h^4$')

    ax.grid(ls=':', which="both")
    ax.set_xlabel(r"$h/\sigma$", fontsize=ftSz2)
    ax.set_ylabel(r"Global error at $ct/L \,=\, {:.3f}$".format(t_wanted), fontsize=ftSz2)
    ax.legend(fontsize=ftSz3)
    
    # fig.savefig("./figures/order.svg", format="svg", bbox_inches='tight')
    plt.show()


def animation_soluce():
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

    with open(filename, 'r') as file:
        scheme = file.readline().strip()
        c, sigma, U_max, L, dt, a = [float(x) for x in file.readline().split()]
        matrix = np.loadtxt(file)

    t = matrix[:, 0]
    u = matrix[:, 1:]
    M, N = u.shape
    h = L / N
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
    _ = FuncAnimation(fig, animate, M, interval=50, blit=True, init_func=init, repeat_delay=3000)

    # to get only one frame at t = i
    # i = M-1 ; init() ; animate(i)
    plt.show()


if __name__ == "__main__":
    
    # plt.rcParams["text.usetex"] = True

    # animation_soluce()

    # plot_soluce()
    plot_wavepacket()
    # plot_soluce_nonuniform()
    # plot_diagnostic()
    # order_convergence()
