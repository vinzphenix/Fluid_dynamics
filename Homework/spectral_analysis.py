import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft, fftshift
from scipy.linalg import eig
from numpy import pi, sqrt, exp, cos

ftSz1, ftSz2, ftSz3 = 20, 15, 12
plt.rcParams["text.usetex"] = False
plt.rcParams['font.family'] = 'monospace'


def fourier():
    U = 1.
    L = 1.
    N = 128
    h = L / N
    sigma = np.array([L / 4., L / 16.])

    j = np.linspace(-N//2, (N-1)//2, N)
    x = np.linspace(-L / 2, L / 2 - h, N)
    k = 2 * pi * j / L

    f_list = U * np.array([exp(-(x / s) ** 2) for s in sigma])
    u_analytic = U * sqrt(np.pi) * np.array([s * exp(-(k * s) ** 2 / 4) for s in sigma])
    u_fft = np.array([np.abs(fftshift(fft(f))) / N for f in f_list])

    fig, axs = plt.subplots(2, 1, figsize=(10, 6), constrained_layout=True, sharex="all")

    # for u_analytic, u_fft, ratio in zip(u_analytic_list, u_fft_list, [4., 16.]):
    axs[0].plot(j, u_analytic[0], "-", lw=3, label='Analytical', alpha=0.85)
    axs[0].plot(j, u_fft[0], ls="", marker="x", label='Numerical')
    axs[1].plot(j, u_analytic[1], "-", lw=3, label='Analytical', alpha=0.85)
    axs[1].plot(j, u_fft[1], ls="", marker="x", label='Numerical')
    
    for ax, ratio in zip(axs, [4., 16.]):
        ax.set_ylim([1.5e-20, 5.])
        ax.set_yscale("log")
        ax.legend(fontsize=ftSz2)
        ax.grid(ls=':')
        ax.set_ylabel(r"$L/\sigma = {:.0f}$".format(ratio), fontsize=ftSz2)

    axs[1].set_xlabel("j", fontsize=ftSz2)
    
    plt.show()


def fourier_packet():
    U, L = 1., 1.
    N = 128
    h = L / N
    s = L / 16.

    j = np.linspace(-N//2, (N-1)//2, N)
    x = np.linspace(-L / 2, L / 2 - h, N)
    k = 2 * pi * j / L
    kp = 2 * pi / N * 16

    f = U * cos(kp * x) * exp(-(x / s) ** 2)
    u_analytic = U / 2. * sqrt(np.pi) * s * exp(-((k - kp) * s) ** 2 / 4)
    u_analytic += U / 2. * sqrt(np.pi) * s * exp(-((k + kp) * s) ** 2 / 4)
    u_fft = np.abs(fftshift(fft(f))) / N

    fig, ax = plt.subplots(1, 1, figsize=(10, 6), constrained_layout=True)

    # for u_analytic, u_fft, ratio in zip(u_analytic_list, u_fft_list, [4., 16.]):
    ax.plot(j, u_analytic, "-", lw=3, label='Analytical', alpha=0.85)
    ax.plot(j, u_fft, ls="", marker="x", label='Numerical')
    
    ax.set_ylim([1.5e-20, 5.])
    ax.set_yscale("log")
    ax.legend(fontsize=ftSz2)
    ax.grid(ls=':')
    ax.set_ylabel(r"fft")
    ax.set_xlabel("j", fontsize=ftSz2)
    
    plt.show()

def stability():
    N = 10
    A = np.zeros((N, N))
    A += np.diag(np.ones(N - 1), 1)
    A -= np.diag(np.ones(N - 1), -1)
    A[0, -1] = -1
    A[-1, 0] = 1
    A *= -1 / 2

    print(A)

    j = np.linspace(0, N - 1, N)
    lambda_exact = -1j * np.sin(2 * pi * j / N)
    lambda_num = eig(A, left=False, right=False)

    print(lambda_exact)
    print(lambda_num)

    fig, ax = plt.subplots(1, 1, figsize=(10, 6), constrained_layout=True)
    ax.plot(np.real(lambda_num), np.imag(lambda_num), 'o', label='Numerical')
    ax.plot(np.real(lambda_exact), np.imag(lambda_exact), '.', label='Analytic')

    ax.set_aspect('equal', 'datalim')
    ax.grid(ls=':')
    ax.legend()
    plt.show()
    return


if __name__ == "__main__":
    # fourier()
    fourier_packet()
    # stability()
