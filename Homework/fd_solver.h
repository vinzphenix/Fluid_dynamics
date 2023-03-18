#ifndef fd_solver_h
#define fd_solver_h

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define SAVE 2  // 0:no save, 1:general for anim, 2:specific file for plots

// Simulation parameters
#define SCHEME_A 'E'
#define SCHEME_B '2'

#define N     128        // nb of points
#define TEND  1.         // final time [s]
#define L     1.         // length [m]
#define C     1.         // wave speed [m/s]
#define SIGMA (L / 16.)  // std deviation
#define UMAX  1.         // height of the gaussian function

// #define KP         2. * M_PI / (L / 8.)    // enable wavepacket mode
#define KP         0.                     // disable wavepacket mode
#define A          0.5                    // parameter of the mapping x(xi) : 0. <= xi < 1.
#define CFL        1.                     // Courant–Friedrichs–Lewy condition
// CFL = (E2: 2.828) (E4: 2.061) (E6: 1.783)
//       (I4: 1.632) (I6: 1.421)

char *path = "./data/";
char filename[50];

// RK4C parameters
double ALPHA[4] = {0, 0.5, 0.5, 1.};
double GAMMA[4] = {1. / 6., 1. / 3., 1. / 3., 1. / 6.};

typedef struct data_Sim_alias {
    int M;
    double h, dt;
    double *u, *us, *ul, *du;  // vectors for the RK4C
    double *dg;                // derivative of the mapping at xi_i
    double *x1, *at, *q;       // vectors for the tri-diagonal system
} data_Sim;

void f_eval(data_Sim *sim);

#endif
