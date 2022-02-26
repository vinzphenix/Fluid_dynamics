#ifndef fd_solver_h
#define fd_solver_h

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define SAVE = 1

// char *filename = "./data/solution_E2_32.txt";
char *filename = "./data/solution.txt";

double ALPHA[4] = {0, 0.5, 0.5, 1.};
double GAMMA[4] = {0.16666666666, 0.33333333333, 0.33333333333, 0.16666666666};

typedef struct data_Sim_alias {
    int N, M;
    double c, sigma, U_max, L, h, dt, CFL;
    double *u, *us, *ul, *du, *q;
    char *scheme;
    void (*f_eval)(struct data_Sim_alias *);
} data_Sim;

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

    // faster than modulo, but not as nice
    du[0] = coef * (2 * (u[1] - u[N-1]) - (u[2] - u[N-2]) / 4.);
    du[1] = coef * (2 * (u[2] - u[0]) - (u[3] - u[N-1]) / 4.);
    du[N - 2] = coef * (2 * (u[N-1] - u[N-3]) - (u[0] - u[N-4]) / 4.);
    du[N - 1] = coef * (2 * (u[0] - u[N-2]) - (u[1] - u[N-3]) / 4.);

    for (int i = 2; i < N - 2; i++) {
        du[i] = coef * (2 * (u[i+1] - u[i-1]) - (u[i+2] - u[i-2]) / 4.);
    }
}

void f_E6(data_Sim *sim) {
    double *u = sim->ul;
    double *du = sim->du;
    int N = sim->N;
    double coef = -sim->c / (20 * sim->h);

    // faster than modulo, but not as nice
    int i;
    for (i = 0; i < 3; i++) {
        du[i] = coef * (15 * (u[i+1] - u[(i-1) % N]) - 3 * (u[i+2] - u[(i-2) % N]) + (u[i+3] - u[(i-3) % N]) / 3.);
    }
    for (i = N - 3; i < N; i++) {
        du[i] = coef * (15 * (u[(i+1) % N] - u[i-1]) - 3 * (u[(i+2) % N] - u[i-2]) + (u[(i+3) % N] - u[i-3]) / 3.);
    }

    for (i = 3; i < N - 3; i++) {
        du[i] = coef * (15 * (u[i+1] - u[i-1]) - 3 * (u[i+2] - u[i-2]) + (u[i+3] - u[i-3]) / 3.);
    }
}

void f_I4(data_Sim *sim) {
    int N = sim->N;
    double *u = sim->ul;
    double *du = sim->du;
    double *q = sim->q;
    double coef = -3 * sim->c / (4 * sim->h);

    q[0] = coef * (u[1] - u[N-1]);
    q[N - 1] = coef * (u[0] - u[N-2]);

    for (int i = 1; i < N - 1; i++) {
        q[i] = coef * (u[i+1] - u[i-1]);
    }

    // du is the solution of the tridiagonal system after this call
    solve_period_3diag(N, 1., 0.25, 0.25, du, q);
}

void f_I6(data_Sim *sim) {
    int N = sim->N;
    double *u = sim->ul;
    double *du = sim->du;
    double *q = sim->q;
    double coef = -1 * sim->c / (9 * sim->h);

    q[0] = coef * (7 * (u[1] - u[N-1]) + (u[2] - u[N-2]) / 4.);
    q[1] = coef * (7 * (u[2] - u[0]) + (u[3] - u[N-1]) / 4.);
    q[N-2] = coef * (7 * (u[N-1] - u[N-3]) + (u[0] - u[N-4]) / 4.);
    q[N-1] = coef * (7 * (u[0] - u[N-2]) + (u[1] - u[N-3]) / 4.);
    
    for (int i = 1; i < N - 1; i++) {
        q[i] = coef * (7 * (u[i+1] - u[i-1]) + (u[i+2] - u[i-2]) / 4.);
    }

    // du is the solution of the tridiagonal system after this call
    solve_period_3diag(N, 1., 1./3., 1./3., du, q);
}

#endif