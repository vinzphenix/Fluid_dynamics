import numpy as np 
import matplotlib.pyplot as plt
import scipy.sparse as sp

ftSz1, ftSz2, ftSz3 = 20, 17, 15
plt.rcParams['font.family'] = 'monospace'

dirPath = "../data"

def plot_matrix():
    matrix = np.loadtxt("matrix.txt")
    row = matrix[:, 0].astype(int)
    col = matrix[:, 1].astype(int)
    val = matrix[:, 2]

    A = sp.coo_matrix((val, (row, col)))

    fig, ax = plt.subplots(1, 1, figsize=(8, 6), constrained_layout=True)
    ax.spy(A)

    plt.show()


def read_data(filename):
    with open(filename, "r") as f:
        n = int(f.readline())
        values = np.loadtxt(f)

    return n, values.T

if __name__ == "__main__":

    n, w = read_data(f"{dirPath}/test_poisson.txt")

    L, H = 15, 5
    h = 1. / n

    # middle of cell : p value of MAC mesh
    xm = np.linspace(h/2., L-h/2., L * n)
    ym = np.linspace(h/2., H-h/2., H * n)
    
    # corners of cell : w value of MAC mesh
    x = np.linspace(0, L, L*n + 1)
    y = np.linspace(0, H, H*n + 1)

    xx_, yy_ = np.meshgrid(x, y)
    xx, yy = np.meshgrid(xm, ym)

    # want u and v at corners of boxes
    # want w and p at center of boxes
    # w_full = np.sin(np.pi * xx_ / L) * np.sin(2 * np.pi * yy_ / H)

    # w = (w_full[:-1, :-1] + w_full[1:, :-1] + w_full[:-1, 1:] + w_full[1:, 1:]) / 4.
    # u = +1./H * np.sin(np.pi * xx_ / L) ** 2 * np.sin(2 * np.pi * yy_ / H)
    # v = -1./L * np.sin(2 * np.pi * xx_ / L) * np.sin(np.pi * yy_ / H) ** 2
    
    mask = (3 < xx) * (xx < 8) * (2 < yy) * (yy < 3)
    # w[mask] = 0.5

    fig = plt.figure(figsize=(11, 9), constrained_layout=True)
    ax = fig.add_subplot(projection='3d')
    
    w = np.ma.masked_array(w, mask)

    vorticity = ax.plot_surface(xx, yy, w, cmap="Spectral")
    ax.set_xlabel("x")
    ax.set_ylabel("y")

    plt.show()


