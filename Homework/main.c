#include "thomas.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

extern double dnrm2_(int *N, double *X, int *INCX);                                   // Norme 2
extern double ddot_(int *N, double *DX, int *INCX, double *DY, int *INCY);            // x dot y
extern void daxpy_(int *N, double *DA, double *DX, int *INCX, double *DY, int *INCY); // y = a*x + y
extern void dcopy_(int *N, double *DX, int *INCX, double *DY, int *INCY);             // copy x to y

typedef struct {
    int N;
    int M;
    double c;
    double sigma;
    double U_max;
    double L;
    double h;
    double dt;
    double CFL;
    double *u;
    double *us;
    double *ul;
    double *du;
    double alpha[4];
    double gamma[4];
} data_Sim;

void init_data_sim(data_Sim *sim, int N, int M, double c, double sigma, double U_max, double L) {
    sim->N = N;
    sim->M = M;
    sim->c = c;
    sim->sigma = sigma;
    sim->U_max = U_max;
    sim->L = L;
    sim->h = L / N;
    sim->CFL = 1.5;
    sim->dt = sim->CFL * sim->h / c;

    sim->u = (double *)malloc(4 * N * sizeof(double));
    sim->ul = sim->u + N;
    sim->us = sim->u + 2 * N;
    sim->du = sim->u + 3 * N;

    double alpha[4] = {0, 0.5, 0.5, 1.};
    double gamma[4] = {0.16666666666, 0.33333333333, 0.33333333333, 0.16666666666};
}

void free_data_sim(data_Sim *sim) {
    free(sim->u);
    free(sim);
}

void set_u_gaussian(data_Sim *sim) {
    double x;
    for (int i = 0; i < sim->N; i++) {
        x = -sim->L / 2 + i * sim->h;
        sim->u[i] = sim->U_max * exp(-x * x / (sim->sigma * sim->sigma));
    }
}

void f_E2(data_Sim *sim) {
    double *u = sim->ul;
    double *du = sim->du;
    int N = sim->N;
    double coef = -sim->c / (2 * sim->h);

    du[0] = coef * (u[1] - u[N - 1]);
    du[N - 1] = coef * (u[0] - u[N - 2]);

    for (int i = 1; i < N - 1; i++) {
        du[i] = coef * (u[i + 1] - u[i - 1]);
    }
}

void f_E4(data_Sim *sim) {
    double *u = sim->ul;
    double *du = sim->du;
    int N = sim->N;
    double coef = -sim->c / (3 * sim->h);

    du[0] = coef * (2 * (u[1] - u[N - 1]) - (u[2] - u[N - 2]) / 4.);
    du[1] = coef * (2 * (u[2] - u[0]) - (u[3] - u[N - 1]) / 4.);
    du[N - 2] = coef * (2 * (u[N - 1] - u[N - 3]) - (u[0] - u[N - 4]) / 4.);
    du[N - 1] = coef * (2 * (u[0] - u[N - 2]) - (u[1] - u[N - 3]) / 4.);

    for (int i = 2; i < N - 2; i++) {
        du[i] = coef * (2 * (u[i + 1] - u[i - 1]) - (u[i + 2] - u[i - 2]) / 4.);
    }
}

void f_I4(data_Sim *sim) {
    double *u = sim->u;
    double *ul = sim->ul;
    double *du = sim->du;
    int N = sim->N;
    double coef = -3 * sim->c / (4 * sim->h);

    ul[0] = coef * (u[1] - u[N - 1]);
    ul[N - 1] = coef * (u[0] - u[N - 2]);
    for (int i = 0; i < N; i++) {
        ul[i] = coef * (u[i + 1] - u[i - 1]);
    }

    // du contains the right values after this call
    solve_period_3diag(N, 1., 0.25, 0.25, du, ul);
}

void save_array(data_Sim *sim, int t) {
    FILE *ptr;
    if (t == 0) ptr = fopen("./data/solution", "w");
    else ptr = fopen("./data/solution", "a");
    
    fprintf(ptr, "%.5lf ", t * sim->dt);
    for (int i = 0; i < sim->N; i++) {
        fprintf(ptr, "%.5e ", sim->u[i]);
    }

    fprintf(ptr, "\n");
    fclose(ptr);
}

void RK4C(data_Sim *sim, void (*f)(data_Sim *)) {
    int N = sim->N;
    double *u = sim->u;
    double *us = sim->us;
    double *ul = sim->ul;
    double *du = sim->du;

    int i, k, t;

    for (t = 0; t < sim->M; t++) {
    
        memset(us, u, N);
        for (k = 0; k < 4; k++) {
            for (i = 0; i < N; i++) {
                ul[i] = us[i] + sim->alpha[k] * sim->dt * du[i];
            }
            f(sim);
            for (i = 0; i < N; i++) {
                u[i] += sim->gamma[k] * sim->dt * du[i];
            }
        }

        save_array(sim, t);
    }

    return;
}

int main(int argc, char *argv[]) {
    data_Sim *simulation = (data_Sim *)malloc(sizeof(data_Sim));
    init_data_sim(simulation, 100, 500, 1., 1. / 16., 1., 1.);
    printf("init done\n");

    set_u_gaussian(simulation);
    printf("set gaussian done\n");

    void (*f_scheme)(data_Sim *);
    if (1) {
        f_scheme = &f_E2;
    } else if (0){
        f_scheme = &f_E4;
    } else {
        f_scheme = &f_I4;
    }

    free_data_sim(simulation);
    return 0;
}