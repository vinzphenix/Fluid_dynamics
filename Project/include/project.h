#ifndef _PROJECT_H_
#define _PROJECT_H_

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>


// Main dimensionless numbers
#define RE 500.            // Reynolds number of the simulation
#define PR 0.7             // Prandtl = nu / alpha
#define GR 100000.         // Grashof = beta (T1-T0) g L^3 / nu^2  // 1000000
#define EC 0.              // Eckert


// Oscillation parameters
#define ALPHA 0.5           // Amplitude of the horizontal oscillation VELOCITY
#define STROUHAL (1. / 3.)  // Frequency of the horizontal oscillation
#define SIWNG_START 0.      // Starting time of the horizontal oscillation

#define KAPPA_Y 0.02387     // Amplitude of the vertical perturbation POSITION  0.02387  0.15
#define STROUHAL_Y (1./3.)  // Frequency of the vertical perturbation
#define N_CYCLES 1          // Duration of the perturbation
#define SMOOTH 0.           // Delay to reach 63% of the desired amplitude of the vertical oscillation (0 to disable)
#define PERT_START 0.       // Starting time of the perturbation


// Temperature parameters
#define NO_SLIP 0           // Walls at y = 0 and y = H
#define TEMP_MODE 0         // Thermal mode 0: disabled, 1: enabled
#define TMIN -1.                 // Min temperature (upper external wall, upper wall of box)
#define TMAX 1.                  // Max temperature (lower external wall, lower wall of box, lateral sides of box)

#if TEMP_MODE
#define WALL_DIRICHLET 1         // external walls dirichlet (0: no flux, 1: dirichlet)
#define BOX_LFT_RGT_DIRICHLET 0  // left and right sides of box dirichlet
#define BOX_BOT_TOP_DIRICHLET 0  // top and bottom sides of box dirichlet
#endif

// Code parameters
#define ADAPTATIVE_DT 1     // 0: classic increment, 1: check reh and rew to adapt time step (but not that functionnal)
#define USE_ADI 0           // 0: classic scheme, 1: solve using ADI method  // *boundary conditions ?
#define CONVECTION_MODE 2   // 0: advective form, 1: divergence form, 2: average of both
#define START_AVG 20.


// Box measurements
#define L_ 15
#define H_ 5
#define LBOX 5
#define D_IN 3
#define D_BOT 2


// Stability settings
#define CFL 0.70
#define FOURIER 0.20
#define U0V0 4.

// Misc
#define TEST_POISSON 0
#define FMT "%.5le\n"
#define ABORT_MSG(s) printf("%s\n%s\n", s, "Program aborted.")

// Progress bar
#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 50


typedef struct {
    int nt, nx, ny, n, first_iteration;
    int size_u, size_v, size_p, size_T;
    int i_start[12], i_final[12], j_start[12], j_final[12];
    double h, dt, tsim, dt_stable;
    double tnow, elapsed, save_freq, start_avg;
    double uMesh, vMesh, reh, rew;
    
    double *u_data, *v_data, *p_data, *T_data;
    double *u_avg, *v_avg;
    double **U, **US, **HX, **HX_;
    double **V, **VS, **HY, **HY_;
    double **P, **PHI;
    double **T, **HT, **HT_;
} Sim_data;


void init_Sim_data(Sim_data *sim, int n_input, double dt_input, double tend_input, double save_freq, const char *myPath);
void init_fields(Sim_data *sim);
void save_fields(Sim_data *sim);
void save_diagnostics(Sim_data *sim, int saving);
void free_Sim_data(Sim_data *sim);
void check_boundary(Sim_data *sim);

void set_bd_conditions(Sim_data *sim);
void set_ghost_points(Sim_data *sim);
void set_boundary_temperature(Sim_data *sim);

void compute_convection(Sim_data *sim);
void compute_convection_temperature(Sim_data *sim);

void predictor_step(Sim_data *sim);

void corrector_step(Sim_data *sim);
void corrector_step_temperature(Sim_data *sim);

void swap_next_previous(Sim_data *sim);
void set_mesh_velocity(Sim_data *sim, double t_now);

#endif
