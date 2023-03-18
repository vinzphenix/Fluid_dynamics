import numpy as np
import matplotlib.pyplot as plt

ftSz1, ftSz2, ftSz3 = 20, 15, 12
plt.rcParams['font.family'] = 'monospace'


def stability_region(h):
    sigma = np.linspace(-4.5, 1.5, 400)
    omega = np.linspace(-4.05, 4.05, 600)
    col = plt.get_cmap('jet')
    x, y = np.meshgrid(sigma, omega)
    z = x + 1j * y
    z2 = z * z
    z4 = z2 * z2

    fig, axs = plt.subplots(2, 3, figsize=(12, 8), constrained_layout=True, sharex='all', sharey='all')
    fig.canvas.manager.set_window_title('Taylor of different orders')

    G = 1. + z
    axs[0, 0].contourf(x, y, np.abs(G), np.arange(0., 1. + h, h), cmap=col)
    axs[0, 0].contour(x, y, np.abs(G), np.arange(0., 1. + h, h), colors='black', linewidths=1)
    axs[0, 0].set_title('order 1')

    G += z2 / 2.
    axs[0, 1].contourf(x, y, np.abs(G), np.arange(0., 1. + h, h), cmap=col)
    axs[0, 1].contour(x, y, np.abs(G), np.arange(0., 1. + h, h), colors='black', linewidths=1)
    axs[0, 1].set_title('order 2')

    G += z2 * z / 6.
    axs[0, 2].contourf(x, y, np.abs(G), np.arange(0., 1. + h, h), cmap=col)
    axs[0, 2].contour(x, y, np.abs(G), np.arange(0., 1. + h, h), colors='black', linewidths=1)
    axs[0, 2].set_title('order 3')

    G += z4 / 24.
    axs[1, 0].contourf(x, y, np.abs(G), np.arange(0., 1. + h, h), cmap=col)
    axs[1, 0].contour(x, y, np.abs(G), np.arange(0., 1. + h, h), colors='black', linewidths=1)
    axs[1, 0].set_title('order 4')

    G += z4 * z / 120. + z4 * z2 / 720
    axs[1, 1].contourf(x, y, np.abs(G), np.arange(0., 1. + h, h), cmap=col)
    axs[1, 1].contour(x, y, np.abs(G), np.arange(0., 1. + h, h), colors='black', linewidths=1)
    axs[1, 1].set_title('order 6')

    G += z4 * z2 * z / (720. * 7.) + z4 * z4 / (720. * 7. * 8.)
    c = axs[1, 2].contourf(x, y, np.abs(G), np.arange(0, 1 + h, h), cmap=col)
    axs[1, 2].contour(x, y, np.abs(G), np.arange(0, 1 + h, h), colors='black', linewidths=1)
    axs[1, 2].set_title('order 8')

    cbar = fig.colorbar(c, ax=axs.ravel().tolist())
    cbar.ax.set_ylabel(r'Module of amplification factor G')
    for ax in axs.flatten():
        ax.grid(ls=':')
        ax.axvline(0., color='grey', lw=3, alpha=0.5, zorder=0)
        ax.axhline(0., color='grey', lw=3, alpha=0.5, zorder=0)
        ax.set_xlim(-5.25, 2.25)
        ax.set_ylim(-4.05, 4.05)
        ax.set_aspect('equal')
    plt.show()
    return


def compare_RK43():
    sigma = np.linspace(-9, 1, 2000)
    omega = np.linspace(-7, 7, 500)
    x, y = np.meshgrid(sigma, omega)
    z = x + 1j * y

    fig, axs = plt.subplots(1, 2, figsize=(12, 8), constrained_layout=True)
    fig.canvas.manager.set_window_title('RK comparison')

    func_rk33 = lambda v: 1 + v + np.power(v, 2) / 2. + np.power(v, 3) / 6.
    # func_rk43 = lambda v: 1 + v + np.power(v, 2) / 2. + np.power(v, 3) / 6. + np.power(v, 4) / 18
    func_rk43 = lambda v: 1 + v + np.power(v, 2) / 2. + np.power(v, 3) / 6. + np.power(v, 4) / (24. * 2.4)
    func_rk44 = lambda v: 1 + v + np.power(v, 2) / 2. + np.power(v, 3) / 6. + np.power(v, 4) / 24.

    G_RK33_re = (8 * func_rk33(z / 2) ** 2 - func_rk33(z)) / 7
    G_RK43_re = (8 * func_rk43(z / 2) ** 2 - func_rk43(z)) / 7
    G_RK44_re = (16 * func_rk44(z / 2) ** 2 - func_rk44(z)) / 15

    list_ = [func_rk33(z), func_rk43(z), func_rk44(z)]
    list_re = [G_RK33_re, G_RK43_re, G_RK44_re]
    list_ct = [None, None, None]

    # G_RK33 = 1 + z + np.power(z, 2) / 2. + np.power(z, 3) / 6.
    # G_RK43 = G_RK33 + np.power(z, 4) / 18.
    # G_RK44 = G_RK33 + np.power(z, 4) / 24.
    for i, rho in enumerate(list_):
        axs[0].contourf(x, y, np.abs(rho), np.array([0., 1.]), colors=f'C{i:d}', alpha=0.1)
        list_ct[i] = axs[0].contour(x, y, np.abs(rho), np.array([0., 1.]), colors=f'C{i:d}')

    for i, rho in enumerate(list_re):
        axs[1].contourf(x, y, np.abs(rho), np.array([0., 1.]), colors=f'C{i:d}', alpha=0.1)
        axs[1].contour(x, y, np.abs(rho), np.array([0., 1.]), colors=f'C{i:d}')

    handles = [c.legend_elements()[0][0] for c in list_ct]
    axs[0].legend(handles, ['RK33', 'RK43', 'RK44'], fontsize=ftSz2, loc='upper left')

    axs[0].set_title("Classic Runge-Kutta", fontsize=ftSz1)
    axs[1].set_title("With Richardson extrapolation", fontsize=ftSz1)
    for ax in axs:
        ax.vlines(0, omega[0], omega[-1], color='black', lw=2., alpha=0.5)
        ax.grid(ls=':')
        ax.set_aspect('equal', 'datalim')

    plt.show()
    return


def rk(u, dt, n, f, beta, gamma):  # Runge-kutta scheme for any number of stages
    gamma = np.array(gamma)  # 1d array
    beta = [np.array([])] + [np.array(beta_row) for beta_row in beta]  # list of 1d arrays
    alpha = np.array([np.sum(beta_row) for beta_row in beta])  # 1d array
    n_stages = np.size(gamma)

    for i in range(n):
        ti = i * dt
        k_list = np.empty(n_stages)
        for stage, (beta_row, alpha_k) in enumerate(zip(beta, alpha)):
            u_input = u[i] + np.dot(beta_row, k_list[:stage])
            k_list[stage] = dt * f(u_input, ti + alpha_k * dt)
        u[i + 1] = u[i] + np.dot(gamma, k_list)

    return


def mesure_error():
    fig, ax = plt.subplots(1, 1, figsize=(10, 6), constrained_layout=True, sharex='all', sharey='all')
    fig.canvas.manager.set_window_title('RK error measure')

    test_f = lambda u, t: -u * np.tanh(t)  # Nonlinear
    u_exact = lambda u0, t: u0 / np.cosh(t)
    # test_f = lambda u, t: -u  # Linear
    # u_exact = lambda u0, t: u0 * np.exp(-t)
    T_final, u_init = 2. * np.pi, 2.  # in [-pi, pi]

    N = np.array([20, 50, 100, 200, 500, 1000])
    DT = T_final / N

    phi = -1e-3
    schemes_params = [
        ([[1.]], [1. / 2., 1. / 2.], "RK2 heun"),
        ([[1. / 3.], [0., 2. / 3.]], [1. / 4., 0., 3. / 4.], "RK3 heun"),
        ([[1. / 2.], [0., 1. / 2.], [0., 0., 1.]], [1. / 6., 1. / 3., 1. / 3., 1. / 6.], "RK4C"),
        ([[2. / 3.], [-phi, phi]], [(phi - 1.) / (4. * phi), 3. / 4., 1. / (4 * phi)], "RK3 special 1"),
        ([[2. / 3.], [2. / 3. - phi, phi]], [1. / 4., (3. * phi - 1) / (4. * phi), 1. / (4 * phi)], "RK3 special 2"),
        ([[1. / 2.], [-1. / 6., 2. / 3.], [1. / 3., -1. / 3., 1.]], [1. / 6., 1. / 3., 1. / 3., 1. / 6.], "RK43"),
    ]
    E = np.zeros((len(schemes_params), np.size(N)))

    for i, (beta, gamma, label) in enumerate(schemes_params):
        for j, (n, dt) in enumerate(zip(N, DT)):
            U = np.zeros((n + 1))
            U[0] = u_init
            rk(U, dt, n, test_f, beta, gamma)
            E[i, j] = np.abs(U[-1] - u_exact(U[0], T_final))
            # ax.plot(np.linspace(0., T_final, n + 1), U, label=r'$\Delta t = {{:.3f}}$'.format(T_final / n))
        ax.loglog(DT, E[i], '-o', label=label)

    for i, (Ei, ls) in enumerate(zip(E, ["--", "-.", ":"])):
        order = i + 2
        label = r"$\propto \Delta t^{:d}$".format(order)
        ax.loglog(DT, np.power(DT / DT[-1], order) * Ei[-1], ls=ls, color='k', label=label)

    ax.legend(fontsize=ftSz2)
    ax.set_ylabel(r"$|\:u^h(T) - u(T)\:|$", fontsize=ftSz2)
    ax.set_xlabel(r"$\Delta t$", fontsize=ftSz2)
    ax.grid(ls=':')
    plt.show()


if __name__ == "__main__":
    plt.rcParams["text.usetex"] = False
    # stability_region(0.1)
    # compare_RK43()
    mesure_error()
