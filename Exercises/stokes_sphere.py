import matplotlib.pyplot as plt
import numpy as np
from numpy import cos, sin, pi
from scipy.integrate import solve_ivp

ftSz1, ftSz2, ftSz3 = 20, 15, 12
plt.rcParams["text.usetex"] = False
plt.rcParams['font.family'] = 'monospace'


def dynamics(_, z):
    c_, s_, e_ = cos(z[1]), sin(z[1]), a / z[0]
    ur = U * c_ / 2. * (e_ ** 3 - 3 * e_ + 2)
    uth = -U * s_ / 4. * (-e_ ** 3 - 3 * e_ + 4)
    return np.array([ur, uth / z[0]])


event = lambda _, z: z[0] * cos(z[1]) - x_max
event.terminal = True
event.direction = 1.

if __name__ == "__main__":
    a, U, mu, p0 = 1., 10., 1., 0.
    n_x, n_y = 30, 10

    x_max = 5. * a
    x_1d = np.linspace(-x_max, x_max, n_x)
    y_1d = np.linspace(0, 3. * a, n_y)
    x, y = np.meshgrid(x_1d, y_1d)

    r = np.hypot(x, y)
    th = np.arctan2(y, x)
    eta = a / r
    sth, cth = sin(th), cos(th)

    u_r = U * cth / 2. * (eta ** 3 - 3 * eta + 2)
    u_th = -U * sth / 4. * (-eta ** 3 - 3 * eta + 4)

    th_o = np.linspace(0, pi, 200)
    r_o = np.linspace(1, 1.5 * x_max / a, 200)
    th_o, r_o = np.meshgrid(th_o, r_o)
    p = -3. / 2. * mu * U / a * cos(th_o) * (a / r_o) ** 2 + p0
    drag = 6. * pi * mu * a * U

    mask = r > a
    # u_r[r < a] = 0.
    # u_th[r < a] = 0.

    ux = u_r * cth - u_th * sth
    uy = u_r * sth + u_th * cth

    init_y = 0.75 * y_1d[1:] ** 2 / np.amax(y_1d)
    init_values = np.c_[np.hypot(-x_max, init_y), np.arctan2(init_y, -x_max)]

    width, height = 12, 8.
    fig, ax = plt.subplots(1, 1, figsize=(width, height), constrained_layout=True)

    for init_value in init_values:
        res = solve_ivp(dynamics, [0, 500.], init_value, 'RK45', events=event, rtol=1.e-8, atol=1.e-7, max_step=0.1)
        this_x = res.y[0] * cos(res.y[1])
        this_y = res.y[0] * sin(res.y[1])
        ax.plot(this_x, this_y, color='black', alpha=0.5)
        ax.plot(this_x, -this_y, color='black', alpha=0.5)

    x_, y_ = x[mask], y[mask]
    ux_, uy_ = ux[mask], uy[mask]
    circle = plt.Circle((0, 0), a, facecolor='w', edgecolor='black', lw=2)
    ax.contour(r_o * cos(th_o), r_o * sin(th_o), p, levels=15, cmap=plt.get_cmap('coolwarm'), linewidths=2)
    ax.contour(r_o * cos(th_o), -r_o * sin(th_o), p, levels=15, cmap=plt.get_cmap('coolwarm'), linewidths=2)
    ax.contourf(r_o * cos(th_o), r_o * sin(th_o), p, levels=15, cmap=plt.get_cmap('coolwarm'), alpha=0.5)
    ax.contourf(r_o * cos(th_o), -r_o * sin(th_o), p, levels=15, cmap=plt.get_cmap('coolwarm'), alpha=0.5)
    ax.quiver(x_, y_, ux_, uy_, np.hypot(ux_, uy_), angles='xy', width=0.003, zorder=2, cmap=plt.get_cmap('viridis'))
    ax.quiver(x_, -y_, ux_, -uy_, np.hypot(ux_, uy_), angles='xy', width=0.003, zorder=2, cmap=plt.get_cmap('viridis'))
    ax.add_patch(circle)

    # ax.grid(ls=':')
    k = 3.
    ax.axis([-k * a * width / height, k * a * width / height, -k * a, k * a])
    ax.set_aspect('equal', 'datalim')
    plt.show()
