import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft, fftshift
from scipy.linalg import eig
from numpy import pi, sqrt, exp, cos, sin

ftSz1, ftSz2, ftSz3 = 20, 17, 15
plt.rcParams["text.usetex"] = False
plt.rcParams['font.family'] = 'monospace'


def fourier(save=False):
    U = 1.
    L = 1.
    N_list = [32, 128]
    ratio = [4., 16.]

    fig, axs = plt.subplots(1, 2, figsize=(10, 5), constrained_layout=True)

    for i, (ax, N, ratio) in enumerate(zip(axs, N_list, ratio)):
        h = L / N
        s = L / ratio
        j = np.linspace(-(N // 2), (N - 1) // 2, N)
        x = np.linspace(-L / 2, L / 2 - h, N)
        k = 2 * pi * j / L

        f = U * exp(-(x / s) ** 2)
        u_analytic = U * sqrt(np.pi) * s * exp(-(k * s) ** 2 / 4)
        u_fft = np.abs(fftshift(fft(f))) / N

        ax.plot(j, u_analytic, "-", lw=3, label='Analytical', alpha=0.85)
        ax.plot(j, u_fft, ls="", marker="x", label='Numerical')

        ax.set_ylim([1.5e-20, 5.])
        ax.set_yscale("log")
        ax.grid(ls=':')
        ax.set_title(r"$L \, / \, \sigma \,=\, {:.0f}$".format(ratio), fontsize=ftSz2)
        ax.set_xlabel("Index j", fontsize=ftSz2)

    axs[0].legend(fontsize=ftSz3)
    if save:
        fig.savefig("./figures/fourier.svg", format="svg", bbox_inches='tight')
    else:
        plt.show()


def fourier_packet(save=False):
    U, L = 1., 1.
    N = 128
    h = L / N
    s = L / 16.

    j = np.linspace(-(N // 2), (N - 1) // 2, N)
    x = np.linspace(-L / 2, L / 2 - h, N)
    k = 2 * pi * j / L
    kp = 2 * pi / L * 16

    f = U * cos(kp * x) * exp(-(x / s) ** 2)
    u_analytic = U / 2. * sqrt(np.pi) * s * exp(-((k - kp) * s) ** 2 / 4)
    u_analytic += U / 2. * sqrt(np.pi) * s * exp(-((k + kp) * s) ** 2 / 4)
    u_fft = np.abs(fftshift(fft(f))) / N

    fig, ax = plt.subplots(1, 1, figsize=(10, 5), constrained_layout=True)

    # for u_analytic, u_fft, ratio in zip(u_analytic_list, u_fft_list, [4., 16.]):
    ax.plot(j, u_analytic, "-", lw=3, label='Analytical FT', alpha=0.85)
    ax.plot(j, u_fft, ls="", marker="x", label='Numerical FFT')

    ax.vlines(-16, 1e-20, 1e5, color='k', ls='--', alpha=0.5)
    ax.vlines(16, 1e-20, 1e5, color='k', ls='--', alpha=0.5, label='Cosine peaks')

    ax.set_ylim([1.5e-20, 5.])
    ax.set_yscale("log")
    ax.legend(fontsize=ftSz3)
    ax.grid(ls=':')
    ax.set_ylabel("Amplitude", fontsize=ftSz2)
    ax.set_xlabel("Index j", fontsize=ftSz2)
    if save:
        fig.savefig("./figures/fourier_packet.svg", format="svg", bbox_inches='tight')
    else:
        plt.show()


def dispersion(save=False):
    fig, axs = plt.subplots(1, 2, figsize=(12, 6), constrained_layout=True, sharex='all')
    n_plot = 400

    x = np.linspace(0, pi, n_plot)
    E2 = sin(x)
    E4 = 4. / 3. * sin(x) - 1. / 6. * sin(2 * x)
    E6 = (45. * sin(x) - 9. * sin(2 * x) + sin(3 * x)) / 30.
    I4 = (3. * sin(x)) / (2. + cos(x))
    I6 = (28. * sin(x) + sin(2. * x)) / (18. + 12 * cos(x))
    axs[0].plot(x / pi, np.divide(E2, x, out=np.ones_like(x), where=x != 0), label='E2', color='C0')
    axs[0].plot(x / pi, np.divide(E4, x, out=np.ones_like(x), where=x != 0), label='E4', color='C1')
    axs[0].plot(x / pi, np.divide(E6, x, out=np.ones_like(x), where=x != 0), label='E6', color='C2')
    axs[0].plot(x / pi, np.divide(I4, x, out=np.ones_like(x), where=x != 0), label='I4', color='C1', ls='--')
    axs[0].plot(x / pi, np.divide(I6, x, out=np.ones_like(x), where=x != 0), label='I6', color='C2', ls='--')
    axs[0].plot(x / pi, np.ones_like(x), label='Exact', color='grey')

    E2 = cos(x)
    E4 = (4. * cos(x) - cos(2 * x)) / 3.
    E6 = (15. * cos(x) - 6. * cos(2 * x) + cos(3 * x)) / 10.
    I4 = 3. * (1 + 2. * cos(x)) / ((2. + cos(x)) ** 2)
    I6 = 1. / 3. * ((28 * sin(x) + sin(2 * x)) * sin(x) + (3 + 2 * cos(x)) * (14 * cos(x) + cos(2 * x))) / (
            (3 + 2 * cos(x)) ** 2)
    axs[1].plot(x / pi, E2, label='E2', color='C0', marker='')
    axs[1].plot(x / pi, E4, label='E4', color='C1', marker='')
    axs[1].plot(x / pi, E6, label='E6', color='C2', marker='')
    axs[1].plot(x / pi, I4, label='I4', color='C1', ls='--', marker='')
    axs[1].plot(x / pi, I6, label='I6', color='C2', ls='--', marker='')
    axs[1].plot(x / pi, np.ones_like(x), label='Exact', color='grey')

    axs[0].set_title("Phase velocity", fontsize=ftSz2)
    axs[0].set_xlabel(r"$kh \: / \: \pi$", fontsize=ftSz2)
    axs[0].set_ylabel(r"$c^* / \: c$", fontsize=ftSz2)
    axs[0].set_aspect("equal", "datalim")

    axs[1].set_title("Group velocity", fontsize=ftSz2)
    axs[1].set_xlabel(r"$kh\: / \: \pi$", fontsize=ftSz2)
    axs[1].set_ylabel(r"$c_g^* \: / \: c$", fontsize=ftSz2)
    for ax in axs:
        ax.grid(ls=':')
        ax.legend(fontsize=ftSz3)

    if save:
        fig.savefig("./figures/dispersion.svg", format="svg", bbox_inches='tight')
    else:
        plt.show()


def dispersion_mosaic(save=False):
    fig, axs = plt.subplot_mosaic([['frq', 'grp'], ['pha', 'grp']], figsize=(10, 7), constrained_layout=True,
                                  sharex=True)
    n_plot = 400

    x = np.linspace(0, pi, n_plot)
    E2 = sin(x)
    E4 = 4. / 3. * sin(x) - 1. / 6. * sin(2 * x)
    E6 = (45. * sin(x) - 9. * sin(2 * x) + sin(3 * x)) / 30.
    I4 = (3. * sin(x)) / (2. + cos(x))
    I6 = (28. * sin(x) + sin(2. * x)) / (18. + 12 * cos(x))

    axs['frq'].plot(x / pi, E2 / pi, label='E2', color='C0')
    axs['frq'].plot(x / pi, E4 / pi, label='E4', color='C1')
    axs['frq'].plot(x / pi, E6 / pi, label='E6', color='C2')
    axs['frq'].plot(x / pi, I4 / pi, label='I4', color='C1', ls='--')
    axs['frq'].plot(x / pi, I6 / pi, label='I6', color='C2', ls='--')
    axs['frq'].plot(x / pi, x / pi, label='Exact', color='grey')

    axs['pha'].plot(x / pi, np.divide(E2, x, out=np.ones_like(x), where=x != 0), label='E2', color='C0')
    axs['pha'].plot(x / pi, np.divide(E4, x, out=np.ones_like(x), where=x != 0), label='E4', color='C1')
    axs['pha'].plot(x / pi, np.divide(E6, x, out=np.ones_like(x), where=x != 0), label='E6', color='C2')
    axs['pha'].plot(x / pi, np.divide(I4, x, out=np.ones_like(x), where=x != 0), label='I4', color='C1', ls='--')
    axs['pha'].plot(x / pi, np.divide(I6, x, out=np.ones_like(x), where=x != 0), label='I6', color='C2', ls='--')
    axs['pha'].plot(x / pi, np.ones_like(x), label='Exact', color='grey')

    E2 = cos(x)
    E4 = (4. * cos(x) - cos(2 * x)) / 3.
    E6 = (15. * cos(x) - 6. * cos(2 * x) + cos(3 * x)) / 10.
    I4 = 3. * (1 + 2. * cos(x)) / ((2. + cos(x)) ** 2)
    I6 = 1. / 3. * ((28 * sin(x) + sin(2 * x)) * sin(x) + (3 + 2 * cos(x)) * (14 * cos(x) + cos(2 * x))) / (
            (3 + 2 * cos(x)) ** 2)
    axs['grp'].plot(x / pi, E2, label='E2', color='C0', marker='')
    axs['grp'].plot(x / pi, E4, label='E4', color='C1', marker='')
    axs['grp'].plot(x / pi, E6, label='E6', color='C2', marker='')
    axs['grp'].plot(x / pi, I4, label='I4', color='C1', ls='--', marker='')
    axs['grp'].plot(x / pi, I6, label='I6', color='C2', ls='--', marker='')
    axs['grp'].plot(x / pi, np.ones_like(x), label='Exact', color='grey')

    axs['frq'].set_title("Circular frequency", fontsize=ftSz2)
    axs['frq'].set_ylabel(r"$\omega^* h / (c\,\pi) $", fontsize=ftSz2)

    axs['pha'].set_title("Phase velocity", fontsize=ftSz2)
    axs['pha'].set_xlabel(r"$kh \: / \: \pi$", fontsize=ftSz2)
    axs['pha'].set_ylabel(r"$c^* / \: c$", fontsize=ftSz2)

    axs['grp'].set_title("Group velocity", fontsize=ftSz2)
    axs['grp'].set_xlabel(r"$kh\: / \: \pi$", fontsize=ftSz2)
    axs['grp'].set_ylabel(r"$c_g^* \: / \: c$", fontsize=ftSz2)

    lines_labels = [axs['frq'].get_legend_handles_labels()]
    lines, labels = [sum(lol, []) for lol in zip(*lines_labels)]
    lgd = fig.legend(lines, labels, labelspacing=2.5, bbox_to_anchor=(0.1, -0.99, 0.80, 1.), mode='expand',
           ncol=6, facecolor='wheat', framealpha=0.25, fancybox=True, fontsize=ftSz3)
    for s, ax in axs.items():
        ax.grid(ls=':')

    if save:
        fig.savefig("./figures/dispersion_bis.svg", format="svg", bbox_extra_artists=(lgd,), bbox_inches='tight')
    else:
        plt.show()


def stability():
    c, L, N = 1., 1., 21
    h = L / N
    A = np.zeros((N, N))
    A += np.diag(np.ones(N - 1), 1)
    A -= np.diag(np.ones(N - 1), -1)
    A[0, -1] = -1
    A[-1, 0] = 1
    A *= -c / (2 * h)

    j = np.linspace(-(N // 2), (N - 1) // 2, N)
    lambda_num = eig(A, left=False, right=False)

    # k = 2 * pi * j / L
    kh = 2 * pi * j / N
    lambda_E2 = -1j * c / h * sin(kh)
    lambda_E4 = -1j * c / h * (4. / 3. * sin(kh) - 1. / 6. * sin(2 * kh))
    lambda_E6 = -1j * c / h * (45. * sin(kh) - 9. * sin(2 * kh) + sin(3 * kh)) / 30.
    lambda_I4 = -1j * c / h * (3. * sin(kh)) / (2. + cos(kh))
    lambda_I6 = -1j * c / h * (28. * sin(kh) + sin(2. * kh)) / (18. + 12 * cos(kh))

    print(np.imag(lambda_I6 * h / c))

    sigma_out, sigma_in  = np.linspace(-3.1, 1.1, 300), np.linspace(-0.01, 0.01, 501)
    omega_out, omega_in  = np.linspace(-3.1, 3.1, 300), np.linspace(-0.1, 0.1, 501)
    sigma, omega = np.sort(np.r_[sigma_out, sigma_in]), np.sort(np.r_[omega_out, omega_in])
    s, w = np.meshgrid(sigma, omega)
    z = s + 1j * w
    G = np.abs(1. + z + np.power(z, 2) / 2. + np.power(z, 3) / 6. + np.power(z, 4) / 24.)

    fig, axs = plt.subplots(1, 2, figsize=(10, 6), constrained_layout=True)
    for ax in axs:
        ax.plot(h / c * np.real(lambda_num), h / c * np.imag(lambda_num), 'o', label='Numerical')
        ax.plot(h / c * np.real(lambda_E2), h / c * np.imag(lambda_E2), '.', label='E2')
        ax.plot(h / c * np.real(lambda_E4), h / c * np.imag(lambda_E4), '.', label='E4')
        ax.plot(h / c * np.real(lambda_E6), h / c * np.imag(lambda_E6), '.', label='E6')
        ax.plot(h / c * np.real(lambda_I4), h / c * np.imag(lambda_I4), '.', label='I4')
        ax.plot(h / c * np.real(lambda_I6), h / c * np.imag(lambda_I6), '.', label='I6')

        # ax.contourf(z.real, z.imag, G, np.array([0., 1.]), cmap=plt.get_cmap('jet'))
        ax.contour(z.real, z.imag, G, np.array([0., 1.]), colors='black', linewidths=1)
        ax.grid(ls=':')

    axs[0].set_aspect('equal', 'datalim')
    axs[0].set_title('RK44', fontsize=ftSz2)
    axs[1].set_title('RK44 zoom', fontsize=ftSz2)
    axs[1].set_xlim([-0.009, 0.009])
    axs[1].legend(fontsize=ftSz3)
    plt.show()


if __name__ == "__main__":
    save_global = False
    plt.rcParams["text.usetex"] = save_global

    # print("{:10s}".format("Plot 1 / 3"), end="\r")
    # fourier(save_global)
    # print("{:10s}".format("Plot 2 / 3"), end="\r")
    # fourier_packet(save_global)
    # print("{:10s}".format("Plot 3 / 3"), end="\r")
    # dispersion_mosaic(save_global)
    # print("{:10s}".format("Job done"))

    # dispersion(save_global)
    stability()
