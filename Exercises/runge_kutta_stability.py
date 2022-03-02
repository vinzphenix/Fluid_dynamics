import numpy as np
import matplotlib.pyplot as plt

ftSz1, ftSz2, ftSz3 = 20, 15, 12
plt.rcParams["text.usetex"] = False
plt.rcParams['font.family'] = 'monospace'


def stability_region(h):
    t1 = np.linspace(-5, 3, 1000)
    t2 = np.linspace(-4, 4, 1000)
    col = plt.get_cmap('jet')
    x, y = np.meshgrid(t1, t2)
    z = x + 1j * y
    z2 = z * z
    # gain = array([1 + z, z2 / 2, z2 * z / 6, z2 * z2 / 24, z2 * z2 * z / 120, z2 * z2 * z2 / 720])

    fig, axs = plt.subplots(2, 3, figsize=(12, 8), constrained_layout=True, sharex='all', sharey='all')
    fig.canvas.manager.set_window_title('Taylor of different orders')

    G = abs(1 + z)
    axs[0, 0].contourf(x, y, G, np.arange(0, 1 + h, h), cmap=col)
    axs[0, 0].contour(x, y, G, np.arange(0, 1 + h, h), colors='black', linewidths=1)
    axs[0, 0].set_title('order 1')

    G = abs(1 + z + z2 / 2)
    axs[0, 1].contourf(x, y, G, np.arange(0, 1 + h, h), cmap=col)
    axs[0, 1].contour(x, y, G, np.arange(0, 1 + h, h), colors='black', linewidths=1)
    axs[0, 1].set_title('order 2')

    G = abs(1 + z + z2 / 2 + z2 * z / 6)
    axs[0, 2].contourf(x, y, G, np.arange(0, 1 + h, h), cmap=col)
    axs[0, 2].contour(x, y, G, np.arange(0, 1 + h, h), colors='black', linewidths=1)
    axs[0, 2].set_title('order 3')

    G = abs(1 + z + z2 / 2 + z2 * z / 6 + z2 * z2 / 24)
    axs[1, 0].contourf(x, y, G, np.arange(0, 1 + h, h), cmap=col)
    axs[1, 0].contour(x, y, G, np.arange(0, 1 + h, h), colors='black', linewidths=1)
    axs[1, 0].set_title('order 4')

    # G = abs(1+z+z2/2+z2*z/6+z2*z2/24+z2*z2*z/120)
    G = abs(1 + z + z ** 2 / 2 + z ** 3 / 6 + z ** 4 / 24 + z ** 5 / 120 + z ** 6 / 720 +
            0 * z ** 7 / (720 * 7) + 0 * z2 ** 4 / (720 * 7 * 8))
    axs[1, 1].contourf(x, y, G, np.arange(0, 1 + h, h), cmap=col)
    axs[1, 1].contour(x, y, G, np.arange(0, 1 + h, h), colors='black', linewidths=1)
    axs[1, 1].set_title('order 6')

    G = abs(1 + z + z ** 2 / 2 + z ** 3 / 6 + z ** 4 / 24 + z ** 5 / 120 + z ** 6 / 720 +
            z ** 7 / (720 * 7) + z2 ** 4 / (720 * 7 * 8))
    c = axs[1, 2].contourf(x, y, G, np.arange(0, 1 + h, h), cmap=col)
    axs[1, 2].contour(x, y, G, np.arange(0, 1 + h, h), colors='black', linewidths=1)
    axs[1, 2].set_title('order 8')

    # fig.colorbar(c)
    cbar = fig.colorbar(c, ax=axs.ravel().tolist())
    cbar.ax.set_ylabel(r'module of Amplification factor')
    for ax in axs.flatten():
        ax.grid(ls=':')
        # ax.set_aspect('equal', 'datalim')

    plt.show()


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


def rk44(u, dt, n, A):
    for i in range(n):
        k1 = dt * np.dot(A, u[i])
        k2 = dt * np.dot(A, u[i] + 1 / 2 * k1)
        k3 = dt * np.dot(A, u[i] + 1 / 2 * k2)
        k4 = dt * np.dot(A, u[i] + k3)
        u[i + 1] = u[i] + 1 / 6 * k1 + 1 / 3 * k2 + 1 / 3 * k3 + 1 / 6 * k4
    return


def rk43_bad(u, dt, n, A):
    for i in range(n):
        k1 = dt * np.dot(A, u[i])
        k2 = dt * np.dot(A, u[i] + 1 / 2 * k1)
        k3 = dt * np.dot(A, u[i] - 1 / 6 * k1 + 2 / 3 * k2)
        k4 = dt * np.dot(A, u[i] + 1 / 3 * k1 - 1 / 3 * k2 + k3)
        u[i + 1] = u[i] + 1 / 6 * k1 + 1 / 3 * k2 + 1 / 3 * k3 + 1 / 6 * k4
    return


def rk43(u, dt, n, A):
    for i in range(n):
        k1 = dt * np.dot(A, u[i])
        k2 = dt * np.dot(A, u[i] + 1 / 2 * k1)
        k3 = dt * np.dot(A, u[i] + 1 / 2 * k2)
        k4 = dt * np.dot(A, u[i] + 7 / 12 * k2 + 5 / 12 * k3)
        u[i + 1] = u[i] + 1 / 6 * k1 + 1 / 3 * k2 + 1 / 3 * k3 + 1 / 6 * k4
    return


def rk33(u, dt, n, A):
    for i in range(n):
        k1 = dt * np.dot(A, u[i])
        k2 = dt * np.dot(A, u[i] + 1 / 3 * k1)
        k3 = dt * np.dot(A, u[i] + 2 / 3 * k2)
        u[i + 1] = u[i] + 1 / 4 * k1 + 3 / 4 * k3
    return


def mesure_error():
    fig, ax = plt.subplots(1, 1, figsize=(10, 6), constrained_layout=True, sharex='all', sharey='all')
    fig.canvas.manager.set_window_title('RK error measure')

    # T_final = 1.
    N = np.array([5, 10, 50, 100, 200, 500, 1000])
    # E = []  # np.zeros_like(N)
    E = np.zeros(len(N))

    for i, n in enumerate(N):
        U = np.zeros((n + 1, 2))
        U[0] = np.array([1, 0])
        rk43(U, 1 / n, n, np.array([[0, 1], [-1, 0]]))
        this_error = U[-1] - np.array([np.cos(1), -np.sin(1)])
        E[i] = np.linalg.norm(this_error, ord=np.inf)

    ax.loglog(N, E, '-o')
    ax.loglog(N, np.power(1 / N, 3.), ls='--')
    ax.grid(ls=':')
    plt.show()


def call_rk():
    N = 100
    DT = 0.01
    U = np.zeros((N + 1, 2))
    U[0] = np.array([1, 0])
    U1 = np.copy(U)
    U2 = np.copy(U)
    rk44(U, DT, N, np.array([[0, 10], [-10, 0]]))
    rk43(U1, DT, N, np.array([[0, 10], [-10, 0]]))
    rk33(U2, DT, N, np.array([[0, 10], [-10, 0]]))

    fig, ax = plt.subplots(1, 1, figsize=(10, 6), constrained_layout=True, sharex='all', sharey='all')
    fig.canvas.manager.set_window_title('RK comparison linear')
    ax.plot(np.linspace(0, N * DT, N + 1), U[:, 0], label='rk44')
    ax.plot(np.linspace(0, N * DT, N + 1), U1[:, 0], label='rk43')
    ax.plot(np.linspace(0, N * DT, N + 1), U2[:, 0], label='rk33')
    ax.legend(fontsize=ftSz2)
    plt.show()


if __name__ == "__main__":
    # stability_region(0.1)
    # compare_RK43()
    # call_rk()
    mesure_error()
