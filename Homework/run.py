import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import os

ftSz1, ftSz2, ftSz3 = 20, 15, 12
plt.rcParams["text.usetex"] = False
plt.rcParams['font.family'] = 'monospace'

filename = "./data/solution.txt"


def plot_function():
    def init():
        exact.set_data(x_plot, f(x_plot, 0))
        time_text.set_text(time_template.format(0))
        line.set_data(x, u[0, :])

        return tuple([line, exact, time_text])

    def animate(t_idx):
        exact.set_data(x_plot, f(x_plot, t_idx * dt))
        time_text.set_text(time_template.format(t_idx * dt))
        line.set_data(x, u[t_idx, :])

        return tuple([line, exact, time_text])

    n_plot = 200
    x_plot = np.linspace(-L / 2., L / 2., n_plot)
    x = np.linspace(-L/2., L/2. - h, N)

    fig, ax = plt.subplots(1, 1, figsize=(10, 6), constrained_layout=True)
    ax.grid(ls=':')
    ax.set_xlabel("x", fontsize=ftSz2)
    ax.set_ylabel("u(x,t)", fontsize=ftSz2)

    time_template = r'$t = {:.2f} [s]$'
    time_text = ax.text(0.03, 0.92, '', fontsize=17, transform=ax.transAxes)

    line, = ax.plot([], [], ls='-', marker='.', color='C0', label='Numerical solution')
    exact, = ax.plot(x_plot, f(x_plot, 0), color='C1', alpha=0.5, lw=5, zorder=0, label='Analytic solution')
    ax.legend(fontsize=ftSz2, loc='upper right')

    # to animate
    _ = FuncAnimation(fig, animate, M, interval=50, blit=True, init_func=init, repeat_delay=5000)

    # to get only one frame at t = i
    # i = M-1 ; init() ; animate(i)

    plt.show()


if __name__ == "__main__":

    with open(filename, 'r') as file:
        c, sigma, U_max, L, dt = [float(x) for x in file.readline().split()]
        matrix = np.loadtxt(file)

    t = matrix[:, 0]
    u = matrix[:, 1:]
    M, N = u.shape
    h = L / N

    f = lambda x_, t_: U_max * np.exp(-np.power((np.fmod(np.abs(x_ - c * t_ + L / 2), L) - L / 2) / sigma, 2))

    plot_function()
