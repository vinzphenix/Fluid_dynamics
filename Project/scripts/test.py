import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys, os

from tqdm import tqdm
from scipy.integrate import cumulative_trapezoid
from scipy.interpolate import RegularGridInterpolator
from time import perf_counter_ns
from matplotlib.animation import FuncAnimation, PillowWriter, FFMpegWriter
from matplotlib.patches import Rectangle
from mpl_toolkits.axes_grid1 import make_axes_locatable


# n = 40
# nx = 15 * n
# ny = 5 * n
# h = 1 / n

# u = np.loadtxt("test_u.txt").reshape(nx+1, ny+2)
# v = np.loadtxt("test_v.txt").reshape(nx+2, ny+1)


# L, H, lbox, din, dbot = 15, 5, 5, 3, 2

# ###############################################################
# st = perf_counter_ns()

# psi1 = np.zeros((nx+1, ny+1))

# for i in range(3*n+1):
#     for j in range(ny):
#         psi1[i, j+1] = psi1[i, j] + h * (-u[i, j+1])

# for i in range(3*n+1, 8*n+1):
#     for j in range(2*n):
#         psi1[i, j+1] = psi1[i, j] + h * (-u[i, j+1])
    
#     psi1[i, 3*n] = psi1[i, 2*n]
#     for j in range(3*n, ny):
#         psi1[i, j+1] = psi1[i, j] + h * (-u[i, j+1])

# for i in range(8*n+1, nx+1):
#     for j in range(ny):
#         psi1[i, j+1] = psi1[i, j] + h * (-u[i, j+1])
# time_exec = perf_counter_ns() - st

# print(f"time method 1 : {time_exec * 1e-6:.3f} ms")
# ###############################################################

# st = perf_counter_ns()

# psi = np.zeros((nx+1, ny+1))

# iwl, iwr = din*n+1, (din+lbox)*n+1
# jwb, jwa = dbot*n+1, (dbot+1)*n+1

# psi[:iwl, 1:] = h * np.cumsum(-u[:iwl, 1:-1], axis=1)                 # block 1
# psi[iwl:iwr, 1:jwb] = h * np.cumsum(-u[iwl:iwr, 1:jwb], axis=1)       # block 2
# psi[iwl:iwr, jwa-1] = psi[iwl:iwr, dbot*n]                            # block 3 
# psi[iwl:iwr, jwa:] = psi[iwl:iwr, jwa-1, np.newaxis] + h * np.cumsum(-u[iwl:iwr, jwa:-1], axis=1)
# psi[iwr:, 1:] = h * np.cumsum(-u[iwr:, 1:-1], axis=1)  # block 4

# time_exec = perf_counter_ns() - st

# print(f"time method 2 : {time_exec * 1e-6:.3f} ms")


# # diff = psi[:din*n+1] - psi_fast[:din*n+1]
# # diff = psi[din*n+1:(din+lbox)*n+1, 1:dbot*n+1] - psi_fast[din*n+1:(din+lbox)*n+1, 1:dbot*n+1]
# # diff = psi[din*n+1:(din+lbox)*n+1, (dbot+1)*n+1:] - psi_fast[din*n+1:(din+lbox)*n+1, (dbot+1)*n+1:]
# # diff = psi[(din+lbox)*n+1:] - psi_fast[(din+lbox)*n+1:]


# print(np.amax(np.abs(psi - psi1)))


# ###############################################################



# x, y = np.meshgrid(np.linspace(0, 15, nx+1), np.linspace(0, 5, ny+1))

# mask = (3 <= x) * (x <= 8) * (2 <= y) * (y <= 3)
# psi = np.ma.masked_array(psi.T, mask)

# psi_min = np.amin(psi)
# psi_max = np.amax(psi)
# psi_avg = (psi_max + psi_min) / 2.
# psi_dif = psi_max - psi_min

# scaling, nStreamLines = 0.75, 20
# levels = np.linspace(psi_min + 0.01 * psi_dif, psi_max - 0.01 * psi_dif, nStreamLines)
# levels = levels - scaling * psi_dif / (2. * np.pi) * np.sin(2. * np.pi / psi_dif * (levels - psi_avg))

# fig, ax = plt.subplots(1, 1, figsize=(10, 5))
# fig.tight_layout()
# rect = Rectangle((3, 2), 5, 1, edgecolor="none", facecolor="black", zorder=10)
# ax.add_patch(rect)
# # ct = ax.pcolormesh(x, y, psi.T)
# ct = ax.contour(x, y, psi, levels=levels)
# cbar = fig.colorbar(ct)

# # ax.plot(x.flatten(), y.flatten(), ls="", marker=".", markersize=1)

# plt.show()


# sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), 'run.py')))
print(sys.path)

from run import read_data_init


