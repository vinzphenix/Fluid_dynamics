#ifndef adi_h
#define adi_h

#include "project.h"

typedef struct {
    int m, n, size;
    double *a, *b, *c, *q;
    double **Q_ux, **Q_vx, **Q_Tx;
} ADI_data;


void init_adi_solver(Sim_data *sim, ADI_data *adi);
void solve_thomas(int size, double *a, double *b, double *c, double *SOL);
void predictor_step_u_adi(Sim_data *sim, ADI_data *adi);
void predictor_step_v_adi(Sim_data *sim, ADI_data *adi);
void free_adi_solver(ADI_data *adi);

#endif
