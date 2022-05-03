import numpy as np

path_res = "/media/vincelinux/LaCie/LMECA2660/results/" + "vertical/" 

filename_stats = path_res + "simu_stats_v0" + ".txt"
filename_params = path_res + "simu_params_v0" + ".txt"

filename_stats_new = path_res + "simu_stats" + ".txt"
filename_params_new = path_res + "simu_params" + ".txt"


with open(filename_params, "r") as f:
    infos = f.readlines()
    line1 = infos[0].split(" ")
    line2 = infos[1].split(" ")
    line3 = infos[2].split(" ")
    line4 = infos[3].split(" ")

    nt, nx, ny, n, save_mod = [int(x) for x in line1]
    RE, tsim, dt, h, L, H, lbox, din, dbot = [float(x) for x in line2]
    alpha, st, swing_start, kappa_y, st_y, pert_start, smooth, n_cycles = [float(x) for x in line3]
    temperature, prandtl, grashof, eckert, tmin, tmax = [float(x) for x in line4]

    n_cycles = int(n_cycles)
    temperature = int(temperature)

with open(filename_params_new, "w") as f:
    # self.nx, self.ny, self.n, self.temperature = params1
    # self.RE, self.T, self.dt, self.h, self.save_freq, self.start_avg = params2
    # self.L, self.H, self.lbox, self.din, self.dbot = params3
    # self.alpha, self.st, self.swing_start, self.kappa_y, self.st_y, self.pert_start, self.smooth, self.n_cycles = params4
    # self.prandtl, self.grashof, self.eckert, self.tmin, self.tmax = params5
    save_freq = save_mod * dt
    start_avg = 20.
    f.write(f"{nx:d} {ny:d} {n:d} {temperature:d}\n")
    f.write(f"{RE:f} {tsim:f} {dt:f} {h:f} {save_freq:f} {start_avg:f}\n")
    f.write(f"{L:f} {H:f} {lbox:f} {din:f} {dbot:f}\n")
    f.write(f"{alpha} {st} {swing_start} {kappa_y} {st_y} {pert_start} {smooth} {n_cycles}\n")
    f.write(f"{prandtl} {grashof} {eckert} {tmin} {tmax}\n")


t_simu = np.linspace(0., tsim, nt + 1)
flag_save = np.zeros(nt+1)
flag_save[::save_mod] = 1
matrix = np.loadtxt(filename_stats)
matrix = np.c_[matrix, t_simu, flag_save]

np.savetxt(filename_stats_new, matrix, fmt=["%.6e", "%.6e", "%.6e", "%.6e", "%.6e", "%.6e", "%.6e", "%d"])
