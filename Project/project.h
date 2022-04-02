#ifndef _PROJECT_H_
#define _PROJECT_H_

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define RE 500.
#define CFL 0.5
#define FOURIER 0.25

#define DEBUG 0
#define OUFLOW_STUPID 0

#define CONVECTION_MODE 0  // 0: advective form, 1: divergence form, 2: both

#define N 50
#define SAVE_MODULO 5
#define FMT "%.4le "

/*
 x_ = x / Hbox
 y_ = y / Hbox
 u_ = u / U0
 v_ = v / U0
 t_ = t / (Hbox/U0)
 p_ = p / (rho U0^2)
*/

#define TEND 0.1
#define L_ 15
#define H_ 5
#define LBOX 5
#define D_IN 3
#define D_BOT 2

#define ALPHA 0.5
#define STROUHAL 0.333333333

// access u and v neighbors
#define LL(w, idx, inc) w[idx-inc]
#define RR(w, idx, inc) w[idx+inc]
#define BB(w, idx, inc) w[idx-1]
#define AA(w, idx, inc) w[idx+1]

// access other field neighbors (reference taken as AL)
#define AL(w, idx, inc) w[idx]       // above left
#define AR(w, idx, inc) w[idx+inc]   // above right
#define BL(w, idx, inc) w[idx-1]     // below left
#define BR(w, idx, inc) w[idx+inc-1] // below right

/*// access u for convection of v
#define LA(w, idx, inc) w[idx-inc+1]
#define RA(w, idx, inc) w[idx+1]
#define LB(w, idx, inc) w[idx-inc]
#define RB(w, idx, inc) w[idx]*/

typedef struct {
    int nt, nx, ny, n;
    int size_u, size_v, size_p;
    int i_start[8], i_final[8], j_start[8], j_final[8];
    double h, dt;
    double uMesh, vMesh;
    
    double *u_data, *v_data, *p_data;
    double **U, **US, **HX, **HX_;
    double **V, **VS, **HY, **HY_;
    double **P, **PHI;

    /*double *u, *u_star, *Hx, *Hx_prev;
    double *v, *v_star, *Hy, *Hy_prev;
    double *p, *phi;*/
} data_Sim;

#endif

