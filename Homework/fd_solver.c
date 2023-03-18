#include "fd_solver.h"

#if ((SCHEME_A == 'E') && (SCHEME_B == '2'))  // error = -h^2/6 d^3u/dx^3
void f_eval(data_Sim *sim) {
    double *u = sim->ul;
    double *du = sim->du;
    double coef = -C / (2. * sim->h);

    du[0] = coef * (u[1] - u[N - 1]);
    du[N - 1] = coef * (u[0] - u[N - 2]);

    for (int i = 1; i < N - 1; i++) {
        du[i] = coef * (u[i + 1] - u[i - 1]);
    }
}

#elif ((SCHEME_A == 'E') && (SCHEME_B == '4'))  // error = h^4/30 d^5u/dx^5
void f_eval(data_Sim *sim) {
    double *u = sim->ul;
    double *du = sim->du;
    double coef = -C / (3. * sim->h);

    // faster than modulo, but not as nice
    du[0] = coef * (2. * (u[1] - u[N-1]) - (u[2] - u[N-2]) / 4.);
    du[1] = coef * (2. * (u[2] - u[0]) - (u[3] - u[N-1]) / 4.);
    du[N - 2] = coef * (2. * (u[N-1] - u[N-3]) - (u[0] - u[N-4]) / 4.);
    du[N - 1] = coef * (2. * (u[0] - u[N-2]) - (u[1] - u[N-3]) / 4.);

    for (int i = 2; i < N - 2; i++) {
        du[i] = coef * (2. * (u[i+1] - u[i-1]) - (u[i+2] - u[i-2]) / 4.);
    }
}

#elif ((SCHEME_A == 'E') && (SCHEME_B == '6'))
void f_eval(data_Sim *sim) {
    double *u = sim->ul;
    double *du = sim->du;
    double coef = -C / (4. * sim->h);

    // faster than modulo everywhere, but not as nice
    int i;
    for (i = 0; i < 3; i++) {
        du[i] = coef * (3. * (u[i+1] - u[(i-1) % N]) - 0.6 * (u[i+2] - u[(i-2) % N]) + (u[i+3] - u[(i-3) % N]) / 15.);
    }
    for (i = N - 3; i < N; i++) {
        du[i] = coef * (3. * (u[(i+1) % N] - u[i-1]) - 0.6 * (u[(i+2) % N] - u[i-2]) + (u[(i+3) % N] - u[i-3]) / 15.);
    }

    for (i = 3; i < N - 3; i++) {
        du[i] = coef * (3. * (u[i+1] - u[i-1]) - 0.6 * (u[i+2] - u[i-2]) + (u[i+3] - u[i-3]) / 15.);
    }
}

#elif ((SCHEME_A == 'D') && (SCHEME_B == '3'))  // error = 2h^3 d^4u/dx^4
void f_eval(data_Sim *sim) {
    // numerical dissipation (not present in other SPATIAL schemes)
    // but still, we use RK4, so there is always some dissipation
    double *u = sim->ul;
    double *du = sim->du;
    double coef = -C / (6. * sim->h);

    int i;
    du[0] = coef * (u[N-2] - 6. * u[N-1] + 3. * u[0] + 2. * u[1]);
    du[1] = coef * (u[N-1] - 6. * u[0] + 3. * u[1] + 2. * u[2]);
    du[N - 1] = coef * (u[N-3] - 6. * u[N-2] + 3. * u[N-1] + 2. * u[0]);

    for (i = 2; i < N - 1; i++) {
        du[i] = coef * (u[i-2] - 6. * u[i-1] + 3. * u[i] + 2. * u[i+1]);
    }
}

#elif ((SCHEME_A == 'I') && (SCHEME_B == '4'))
void f_eval(data_Sim *sim) {
    double *u = sim->ul;
    double *du = sim->du;
    double *q = sim->q;
    double coef = -3. * C / (4. * sim->h);

    q[0] = coef * (u[1] - u[N-1]);
    q[N - 1] = coef * (u[0] - u[N-2]);

    for (int i = 1; i < N - 1; i++) {
        q[i] = coef * (u[i+1] - u[i-1]);
    }

    // du is the solution of the tridiagonal system after this call
    solve_period_3diag(N, 1., 0.25, 0.25, du, q, sim->x1, sim->at);
    //solve_period_3diag(N, 1., 0.25, 0.25, du, q);
}

#elif ((SCHEME_A == 'I') && (SCHEME_B == '6'))
void f_eval(data_Sim *sim) {
    double *u = sim->ul;
    double *du = sim->du;
    double *q = sim->q;
    double coef = -C / (9. * sim->h);

    q[0] = coef * (7. * (u[1] - u[N-1]) + (u[2] - u[N-2]) / 4.);
    q[1] = coef * (7. * (u[2] - u[0]) + (u[3] - u[N-1]) / 4.);
    q[N-2] = coef * (7. * (u[N-1] - u[N-3]) + (u[0] - u[N-4]) / 4.);
    q[N-1] = coef * (7. * (u[0] - u[N-2]) + (u[1] - u[N-3]) / 4.);
    
    for (int i = 1; i < N - 1; i++) {
        q[i] = coef * (7. * (u[i+1] - u[i-1]) + (u[i+2] - u[i-2]) / 4.);
    }

    // du is the solution of the tridiagonal system after this call
    solve_period_3diag(N, 1., 1./3., 1./3., du, q, sim->x1, sim->at);
}
#else
#warning "This scheme is not yet implemented"
#endif
