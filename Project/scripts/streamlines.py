import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from numpy import sin, cos, exp, pi
from scipy.interpolate import LinearNDInterpolator, RegularGridInterpolator
from scipy.integrate import solve_ivp

ftSz1, ftSz2, ftSz3 = 20, 17, 15
plt.rcParams['font.family'] = 'monospace'


def init_stream():
    Q.set_UVC(U[0], V[0])
    time_text.set_text(time_template.format(0))
    
    for j, line in enumerate(lines):
        line.set_data([], [])
    return Q, time_text, *lines

def animate_stream(i):
    tau = t[i]
    
    for j, line in enumerate(lines):
        sol = solve_ivp(f_streamline, [0, 10.], np.array([1e-3, -H/2 + (j+1) / (nLines+1) * H]), args=(tau,), events=event_x)
        line.set_data(sol.y[0], sol.y[1])

    Q.set_UVC(U[i], V[i])
    time_text.set_text(time_template.format(t[i]))
    return Q, time_text, *lines


event_x = lambda s, pos, tau: L - np.abs(2 * pos[0] - L)
event_y = lambda s, pos, tau: H - np.abs(2 * pos[1])
event_x.terminal = True
event_y.terminal = True

def f_streamline(s, position, tau):
    x, y = position
    dx, = u_interp([tau, x, y])
    dy, = v_interp([tau, x, y])
    return np.array([dx, dy])


def init():
    
    animated_items = []

    if kwargs["streak"]:
        for j in range(nLines):
            streaklines_x[j] = -np.ones(nParticles)
            streaklines_y[j] = spots_y[j]*np.ones(nParticles)        
        for j, line in enumerate(lines):
            line.set_data([], [])
        animated_items += lines

    if kwargs["quiver"]:
        Q.set_UVC(U[0], V[0])
        animated_items.append(Q)
        
    # if kwargs["stream"]:
    #     animated_items.append(?)

    if kwargs["vorticity"]:
        pc.set_array(omega[0].flatten())
        animated_items.append(pc)


    time_text.set_text(time_template.format(0))
    animated_items.append(time_text)

    return tuple(animated_items)


def animate(i):
    
    animated_items = []
    
    if kwargs["streak"]:

        for j, (line, x_particles, y_particles) in enumerate(zip(lines, streaklines_x, streaklines_y)):

            # add new particle
            if j == 0 and i % 10 == 0:
                print((i//10) % nParticles)
                x_particles[(i // 10) % nParticles] = 0.
                y_particles[(i // 10) % nParticles] = spots_y[j]

            mask = (0. <= x_particles) * (x_particles <= L) * (-H/2. <= y_particles) * (y_particles <= H/2)

            args = np.c_[t[i] * np.ones(nParticles)[mask], x_particles[mask], y_particles[mask]]

            # Euler explicite
            x_particles[mask] += dt * u_interp(args)
            y_particles[mask] += dt * v_interp(args)

            # RK22
            # x_temp, y_temp = -np.ones(nParticles), np.zeros(nParticles)
            # x_temp[mask] = x_particles[mask] + dt / 2. * u_interp(args)
            # y_temp[mask] = y_particles[mask] + dt / 2. * v_interp(args)
            # mask = (0. <= x_temp) * (x_temp <= L) * (-H/2. <= y_temp) * (y_temp <= H/2)
            # args = np.c_[(t[i] + dt / 2.) * np.ones(nParticles)[mask], x_temp[mask], y_temp[mask]]
            # x_particles[mask] += dt * u_interp(args)
            # y_particles[mask] += dt * v_interp(args)

            x_particles[~mask] = np.nan
            y_particles[~mask] = np.nan
            line.set_data(x_particles, y_particles)
            # line.set_data(x_particles[mask], y_particles[mask])
        
        animated_items += lines

    if kwargs["quiver"]:
        Q.set_UVC(U[i], V[i])
        animated_items.append(Q)

    if kwargs["stream"]:
        ax.cla()
        _, dudy = np.gradient(U[i].T, L / nx, H / ny)
        dvdx, _ = np.gradient(V[i].T, L / nx, H / ny)
        stream = ax.streamplot(x, y, U[i].T, V[i].T, color=np.abs(dvdx - dudy), cmap=plt.get_cmap('viridis'), density=1.)
        animated_items.append(stream)

    if kwargs["vorticity"]:
        pc.set_array(omega[i].flatten())
        animated_items.append(pc)


    time_text.set_text(time_template.format(t[i]))
    animated_items.append(time_text)

    return tuple(animated_items)


if __name__ == "__main__":

    kwargs = {"stream":False, "streak":True, "quiver":True, "vorticity": True}

    nt, nx, ny = 500, 40, 20
    T, L, H = 10, 2., 1.
    x = np.linspace(0., L, nx)
    y = np.linspace(-H/2., H/2., ny)
    t, dt = np.linspace(0., T, nt, retstep=True)
    
    
    X, Y = np.meshgrid(x, y, indexing="ij")

    U = np.einsum("i,j,k -> ijk", 0.1+cos(t)**2, 0.5+0.5*(1 - 2*x/L)**2, 1+cos(pi * y / H))
    V = np.einsum("i,j,k -> ijk", sin(t), 0.5+0.5*(1 - 2*x/L)**2, cos(pi * y / H)**2)

    u_interp = RegularGridInterpolator((t, x, y), U, bounds_error=False, fill_value=0)
    v_interp = RegularGridInterpolator((t, x, y), V, bounds_error=False, fill_value=0)

    nLines, nParticles = 3, nx
    spots_y = np.linspace(-H/2 + H/(nLines+1), H/2 - H/(nLines+1), nLines)
    streaklines_x = [-np.ones(nParticles) for _ in range(nLines)]
    streaklines_y = [spots_y[j]*np.ones(nParticles) for j in range(nLines)]

    fig, ax = plt.subplots(1, 1, figsize=(14, 6), constrained_layout=True)

    if kwargs["quiver"]:
        Q = ax.quiver(X, Y, U[0], V[0], width=0.0015, zorder=2)
    
    if kwargs["streak"]:
        lines = [ax.plot([], [], "-", marker="o", color="yellow", lw=3,markersize=10, alpha=0.75)[0] for i in range(nLines)]

    if kwargs["stream"]:
        1

    if kwargs["vorticity"]:
        dudy = [np.gradient(U[i].T, L / nx, H / ny)[1] for i in range(nt)]
        dvdx = [np.gradient(V[i].T, L / nx, H / ny)[0] for i in range(nt)]
        omega = [dvdx[i] - dudy[i] for i in range(nt)]
        wmin = np.amin(np.array([np.amin(omega[i]) for i in range(nt)]))
        wmax = np.amax(np.array([np.amax(omega[i]) for i in range(nt)]))
        pc = plt.pcolormesh(x, y, omega[0], vmin=wmin, vmax=wmax, cmap=plt.get_cmap("bwr"))  # RdBu

    ax.set_xlim([0., L])
    ax.set_ylim([-0.5 * H, 0.5 * H])
    ax.set_aspect('equal')

    time_template = "t = {:4.2f} [s]"    
    time_text = ax.text(0.8,0.9, "",
                        fontsize=ftSz2,
                        transform=ax.transAxes,
                        bbox=dict(boxstyle="round", fc="wheat", ec="none", alpha=0.85))
    
    _ = FuncAnimation(fig, animate, nt, interval=20, blit=True, init_func=init, repeat_delay=3000)

    # ax.streamplot(x, y, U[idx].T, V[idx].T)

    plt.show()
