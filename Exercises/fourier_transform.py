import numpy as np
import matplotlib.pyplot as plt
from scipy import fft, signal
from matplotlib.animation import FuncAnimation

ftSz1, ftSz2, ftSz3 = 20, 15, 12
plt.rcParams["text.usetex"] = False
plt.rcParams['font.family'] = 'monospace'


def init():
    time_text.set_text(time_template.format(0))
    line1.set_data([], [])
    for line in lines:
        line.set_alpha(alpha)
    return tuple([*lines, line1, time_text])


def animate(i):
    time_text.set_text(time_template.format(i))
    line1.set_data([x_plot], [np.sum(y_fourier[:i+1], axis=0)])
    if i > 0:
        lines[i-1].set_alpha(alpha)
        lines[i-1].set_linewidth(1.5)
    lines[i].set_alpha(1)
    lines[i].set_linewidth(2.5)
    return tuple([*lines, line1, time_text])


if __name__ == "__main__":

    L, nx = 2.5, 33
    n_plot = 1000
    # f1 = lambda x_: 2 / np.pi * np.arcsin(np.sin(6 * np.pi * x_ / L)) + 0.5
    f1 = lambda x_: signal.sawtooth(8 * np.pi * x_ / L)
    f2 = lambda x_: signal.square(4 * np.pi * x_ / L)
    f3 = lambda x_: np.sin(2 * np.pi * x_ / L)
    f4 = lambda x_: np.sin(2 * np.pi * x_ / L) + 0.2 * np.sin(4 * np.pi * x_ / L) + 0.4 * np.sin(6 * np.pi * x_ / L) + 2.
    f5 = lambda x_: np.sin(6 * np.pi * (x_ / L) ** 2.) + np.cos(2 * np.pi * x_ / L)
    f = f1

    x = np.linspace(0, L * (nx - 1) / nx, nx)
    y = f(x)
    A = fft.fft(y)

    x_plot = np.linspace(0, L, n_plot)
    y_fourier = np.zeros((nx // 2 + 1, n_plot))
    y_fourier[0] = A[0].real / nx
    for k in range(1, nx // 2 + 1):
        arg = 2 * np.pi * k * (x_plot / L)
        y_fourier[k] = 2. * (A[k].real * np.cos(arg) - A[k].imag * np.sin(arg)) / nx
        if (nx % 2 == 0) and (k == nx // 2):
            y_fourier[k] /= 2.

    y_sum = np.sum(y_fourier, axis=0)

    fig, axs = plt.subplots(2, 1, figsize=(12, 8), constrained_layout=True)
    axs[0].plot(x, y, 'o', color='C0', label='samples')
    axs[0].plot(x_plot, f(x_plot), '-', color='C0', label='f(x)')
    line1, = axs[0].plot(x_plot, np.zeros_like(x_plot), color='C1', label='fourier series')

    alpha = 0.25
    lines = [axs[1].plot(x_plot, y_fourier[i], alpha=alpha)[0] for i in range(nx // 2 + 1)]

    time_template = r"mode = ${:d}$"
    time_text = axs[0].text(0.05, 0.85, "", transform=axs[0].transAxes, fontsize=ftSz2)

    axs[1].set_xlabel(r"$x$", fontsize=ftSz2)
    axs[0].set_ylabel(r"$f(x)$", fontsize=ftSz2)
    axs[1].set_ylabel(r"Basis functions", fontsize=ftSz2)
    axs[0].grid(ls=':')
    axs[1].grid(ls=':')
    delta = np.amax(y) - np.amin(y)
    axs[0].axis([-L / 25, L, np.amin(y) - delta / 10, np.amax(y) + delta / 10])
    axs[0].legend(loc='lower left', fontsize=ftSz3)

    anim = FuncAnimation(fig, animate, nx // 2 + 1, interval=500, blit=True, init_func=init, repeat_delay=1000)

    plt.show()
