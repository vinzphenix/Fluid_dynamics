import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft, fftshift
from scipy.linalg import eig
from numpy import pi, sqrt, exp

ftSz1, ftSz2, ftSz3 = 20, 15, 12
plt.rcParams["text.usetex"] = False
plt.rcParams['font.family'] = 'monospace'


def fourier():
    U = 1.
    L = 1.
    N = 50
    h = L / N
    s = 1. / 16.

    j = np.linspace(0, N - 1, N)
    x = np.linspace(-L / 2, L / 2 - h, N)
    k = 2 * pi * j / L
    f = U * exp(-(x / s) ** 2)
    u_analytic = U * sqrt(np.pi) * s * exp(-(k * s) ** 2 / 4)

    u_fft = np.abs(fft(f)) / N

    fig, ax = plt.subplots(1, 1, figsize=(10, 6), constrained_layout=True)

    ax.plot(j, 20 * np.log10(u_analytic), "-o", label='analytic')
    ax.plot(j, 20 * np.log10(u_fft), "-o", label='fft')

    ax.legend()
    ax.grid(ls=':')
    ax.set_xlabel("j")
    ax.set_ylabel("dB")
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
    stability()
