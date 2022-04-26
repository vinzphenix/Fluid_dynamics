#ifndef _PROJECT_H_
#define _PROJECT_H_

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>


// #define N_ 40
// #define NT 400
// #define DT 0.002
// #define TSIM 10.           // Final time of the simulation

// Simulation parameters
#define RE 500.            // Reynolds number of the simulation

// Oscillation parameters
#define ALPHA 0.5           // Amplitude of the horizonatal oscillation VELOCITY
#define STROUHAL (1. / 3.)  // Frequency of the horizontal oscillation
#define SIWNG_START 20.    // Starting time of the horizontal oscillation

#define KAPPA_Y 0.05       // Amplitude of the vertical perturbation POSITION
#define STROUHAL_Y (1./3.)  // Frequency of the vertical perturbation
#define PERT_START 20.0      // Starting time of the perturbation
#define N_CYCLES 1          // Duration of the perturbation

// Temperature parameters
#define NO_SLIP 0           // Walls at y = 0 and y = H
#define TEMP_MODE 2         // Thermal mode, 0: no coupling, 1: external walls hot and cold, 2: rectangle heat source
#define PR 2.              // Prandtl = nu / alpha
#define GR 200000.         // Grashof = beta (T1-T0) g L^3 / nu^2  // 1000000
#define NU 0.               // Nusselt
#define EC 0.               // Eckert
#define TMIN -1.            // Min temperature
#define TMAX 0.75             // Max temperature

// Code parameters
#define USE_ADI 0           // 0: classic scheme, 1: solve using ADI method  // *boundary conditions ?
#define CONVECTION_MODE 2   // 0: advective form, 1: divergence form, 2: average of both
#define SAVE 1              // 1 to save, 0 otherwise
#define SAVE_MODULO 50      // save results every ... iteration

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


#define TEST_POISSON 0
#define FMT "%.5le\n"
#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60


typedef struct {
    int nt, nx, ny, n, t;
    int size_u, size_v, size_p, size_T;
    int i_start[12], i_final[12], j_start[12], j_final[12];
    double h, dt, tsim;
    double uMesh, vMesh;
    
    double *u_data, *v_data, *p_data, *T_data;
    double **U, **US, **HX, **HX_;
    double **V, **VS, **HY, **HY_;
    double **P, **PHI;
    double **T, **HT, **HT_;
} Sim_data;


void init_Sim_data(Sim_data *sim, int n_input, double dt_input, double tend_input);
void init_fields(Sim_data *sim);
void save_fields(Sim_data *sim, int t);
void write_diagnostics(Sim_data *sim, int t);
void free_Sim_data(Sim_data *sim);
void check_boundary(Sim_data *sim);

void set_bd_conditions(Sim_data *sim, double **U, double **V);
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
