#include "fd_solver.h"

#if ((SCHEME_A == 'E') && (SCHEME_B == '2'))
void f_eval(data_Sim *sim) {
    double *u = sim->ul;
    double *du = sim->du;
    double coef = -C / (2 * sim->h);

    du[0] = coef * (u[1] - u[N - 1]);
    du[N - 1] = coef * (u[0] - u[N - 2]);

    for (int i = 1; i < N - 1; i++) {
        du[i] = coef * (u[i + 1] - u[i - 1]);
    }
}

#elif ((SCHEME_A == 'E') && (SCHEME_B == '4'))
void f_eval(data_Sim *sim) {
    double *u = sim->ul;
    double *du = sim->du;
    double coef = -C / (3 * sim->h);

    // faster than modulo, but not as nice
    du[0] = coef * (2 * (u[1] - u[N-1]) - (u[2] - u[N-2]) / 4.);
    du[1] = coef * (2 * (u[2] - u[0]) - (u[3] - u[N-1]) / 4.);
    du[N - 2] = coef * (2 * (u[N-1] - u[N-3]) - (u[0] - u[N-4]) / 4.);
    du[N - 1] = coef * (2 * (u[0] - u[N-2]) - (u[1] - u[N-3]) / 4.);

    for (int i = 2; i < N - 2; i++) {
        du[i] = coef * (2 * (u[i+1] - u[i-1]) - (u[i+2] - u[i-2]) / 4.);
    }
}

#elif ((SCHEME_A == 'E') && (SCHEME_B == '6'))
void f_eval(data_Sim *sim) {
    double *u = sim->ul;
    double *du = sim->du;
    double coef = -C / (20 * sim->h);

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

#elif ((SCHEME_A == 'I') && (SCHEME_B == '4'))
void f_eval(data_Sim *sim) {
    double *u = sim->ul;
    double *du = sim->du;
    double *q = sim->q;
    double coef = -3 * C / (4 * sim->h);

    q[0] = coef * (u[1] - u[N-1]);
    q[N - 1] = coef * (u[0] - u[N-2]);

    for (int i = 1; i < N - 1; i++) {
        q[i] = coef * (u[i+1] - u[i-1]);
    }

    // du is the solution of the tridiagonal system after this call
    solve_period_3diag(N, 1., 0.25, 0.25, du, q);
}

#elif ((SCHEME_A == 'I') && (SCHEME_B == '6'))
void f_eval(data_Sim *sim) {
    double *u = sim->ul;
    double *du = sim->du;
    double *q = sim->q;
    double coef = -1 * C / (9 * sim->h);

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

#else
#warning "This scheme is not yet implemented"

#endif
