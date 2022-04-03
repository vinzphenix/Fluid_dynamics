#ifndef adi_h
#define adi_h

#include "project.h"

typedef struct {
    int m, n, size;
    double *a, *b, *c, *q;
    double **A_ux, **B_ux, **C_ux, **Q_ux;
    double **A_uy, **B_uy, **C_uy, **Q_uy;
    double **A_vx, **B_vx, **C_vx, **Q_vx;
    double **A_vy, **B_vy, **C_vy, **Q_vy;
} ADI_data;


void init_adi_solver(Sim_data *sim, ADI_data *adi);
void solve_thomas(int size, double *a, double *b, double *c, double *SOL);
void predictor_step_adi(Sim_data *sim, ADI_data *adi);
void free_adi_solver(ADI_data *adi);

#endif
