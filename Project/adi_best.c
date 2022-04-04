#include "adi.h"
#define MAX(x, y) (((x) > (y)) ? (x) : (y))


void init_adi_solver(Sim_data *sim, ADI_data *adi) {
    int j, size;
    size = MAX(sim->size_u, sim->size_v);
    
    adi->a = (double *)calloc(4 * size, sizeof(double));
    adi->b = adi->a + size;
    adi->c = adi->b + size;
    adi->q = adi->c + size;

    size = (sim->ny + 2) + (sim->ny + 1);
    adi->Q_ux = (double **)malloc(size * sizeof(double *));
    adi->Q_vx = adi->Q_ux + (sim->ny + 2);

    for (j = 0; j < sim->ny+2; j++)
        adi->Q_ux[j] = adi->q + j * (sim->nx+1);
    for (j = 0; j < sim->ny+1; j++)
        adi->Q_vx[j] = adi->q + j * (sim->nx+2);
}


void solve_thomas(int size, double *a, double *b, double *c, double *sol) {
    /*
    Solves the tridiagonal system using the Thomas algorithm. 
    The matrix has the following structure:
     -                 -
     |b c              |
     |a b c            |
     |  a b c          | = M       such that M x = sol
     |        ...      |
     |            a b  |
     -                 -
     
     The solution x is stored in-place, in the vector sol
    */

    int i;
    double w;

    for (i = 1; i < size; i++) {
        w = a[i] / b[i-1];
        b[i] -= w * c[i-1];
        sol[i] -= w * sol[i-1];
    }

    sol[size-1] /= b[size-1];
    for (i = size-2; i >= 0; i--) {
        sol[i] = (sol[i] - c[i] * sol[i+1]) / b[i];
    }
}


void reset_tridiagonal(ADI_data *adi, int size) {
    for (int i = 0; i < size; i++) {
        adi->a[i] = 0.;
        adi->b[i] = 1.;
        adi->c[i] = 0.;
        // adi->q[i] = 0.;
    }
}


void set_system_ux(Sim_data *sim, ADI_data *adi) {
    int i, j, k, idx;
    k = sim->nx + 1;

    // set identity for upper and lower walls
    for (i = 0; i < sim->nx + 1; i++) {
        for (j = 0; j < sim->ny + 2; j += sim->ny+1) {
            idx = j*k + i;
            adi->a[idx] = 0.; adi->b[idx] = 1.; adi->c[idx] = 0.;
            adi->Q_ux[j][i] = sim->U[i][j];
        }
    }

    // set identity for left and right walls
    for (j = 1; j < sim->ny + 1; j++) {
        for (i = 0; i < sim->nx+1; i += sim->nx) {
            idx = j*k + i;
            adi->a[idx] = 0.; adi->b[idx] = 1.; adi->c[idx] = 0.;
            adi->Q_ux[j][i] = sim->U[i][j];
        }
    }

    // set identity for rectangle
    for (i = sim->n * D_IN; i < (D_IN + LBOX) * sim->n + 1; i++) {
        for (j = D_BOT * sim->n + 1; j < (D_BOT + 1) * sim->n + 1; j++) {
            idx = j*k + i;
            adi->a[idx] = 0.; adi->b[idx] = 1.; adi->c[idx] = 0.;
            adi->Q_ux[j][i] = sim->U[i][j];
        }
    }
}

void set_system_uy(Sim_data *sim, ADI_data *adi) {
    int i, j, k, idx;
    k = sim->ny + 2;

    // set identity for upper and lower walls
    for (i = 0; i < sim->nx + 1; i++) {
        for (j = 0; j < sim->ny + 2; j += sim->ny+1) {
            idx = i*k + j;
            adi->a[idx] = 0.; adi->b[idx] = 1.; adi->c[idx] = 0.;
            sim->US[i][j] = sim->U[i][j];
        }
    }

    // set identity for left and right walls
    for (j = 1; j < sim->ny + 1; j++) {
        for (i = 0; i < sim->nx+1; i += sim->nx) {
            idx = i*k + j;
            adi->a[idx] = 0.; adi->b[idx] = 1.; adi->c[idx] = 0.;
            sim->US[i][j] = sim->U[i][j];
        }
    }

    // set identity for rectangle
    for (i = sim->n * D_IN; i < (D_IN + LBOX) * sim->n + 1; i++) {
        for (j = D_BOT * sim->n + 1; j < (D_BOT + 1) * sim->n + 1; j++) {
            idx = i*k + j;
            adi->a[idx] = 0.; adi->b[idx] = 1.; adi->c[idx] = 0.;
            sim->US[i][j] = sim->U[i][j];
        }
    }
}

void set_system_vx(Sim_data *sim, ADI_data *adi) {
    int i, j, k, idx;
    k = sim->nx + 2;

    // set identity for upper and lower walls
    for (i = 0; i < sim->nx + 2; i++) {
        for (j = 0; j < sim->ny + 1; j += sim->ny) {
            idx = j*k + i;
            adi->a[idx] = 0.; adi->b[idx] = 1.; adi->c[idx] = 0.;
            adi->Q_vx[j][i] = sim->U[i][j];
        }
    }

    // set identity for left and right walls
    for (j = 1; j < sim->ny; j++) {
        for (i = 0; i < sim->nx+2; i += sim->nx + 1) {
            idx = j*k + i;
            adi->a[idx] = 0.; adi->b[idx] = 1.; adi->c[idx] = 0.;
            adi->Q_vx[j][i] = sim->U[i][j];
        }
    }

    // set identity for rectangle
    for (i = sim->n * D_IN + 1; i < (D_IN + LBOX) * sim->n + 1; i++) {
        for (j = D_BOT * sim->n; j < (D_BOT + 1) * sim->n + 1; j++) {
            idx = j*k + i;
            adi->a[idx] = 0.; adi->b[idx] = 1.; adi->c[idx] = 0.;
            adi->Q_vx[j][i] = sim->U[i][j];
        }
    }
}

void set_system_vy(Sim_data *sim, ADI_data *adi) {
    int i, j, k, idx;
    k = sim->ny + 1;

    // set identity for upper and lower walls
    for (i = 0; i < sim->nx + 2; i++) {
        for (j = 0; j < sim->ny + 1; j += sim->ny) {
            idx = i*k + j;
            adi->a[idx] = 0.; adi->b[idx] = 1.; adi->c[idx] = 0.;
            sim->VS[i][j] = sim->U[i][j];
        }
    }

    // set identity for left and right walls
    for (j = 1; j < sim->ny; j++) {
        for (i = 0; i < sim->nx+2; i += sim->nx + 1) {
            idx = i*k + j;
            adi->a[idx] = 0.; adi->b[idx] = 1.; adi->c[idx] = 0.;
            sim->VS[i][j] = sim->U[i][j];
        }
    }

    // set identity for rectangle
    for (i = sim->n * D_IN + 1; i < (D_IN + LBOX) * sim->n + 1; i++) {
        for (j = D_BOT * sim->n; j < (D_BOT + 1) * sim->n + 1; j++) {
            idx = i*k + j;
            adi->a[idx] = 0.; adi->b[idx] = 1.; adi->c[idx] = 0.;
            sim->VS[i][j] = sim->U[i][j];
        }
    }
}


void check_adi(Sim_data *sim) {
    int i, j;
    for (i = 0; i < sim->nx+1; i++) {
        for (j = 0; j < sim->ny+2; j++) {
            PRINTF("u[%4d][%4d] = %.5lf  vs us = %.5lf\n", i, j, sim->U[i][j], sim->US[i][j]);
        }
    }
}


void predictor_step_adi(Sim_data *sim, ADI_data *adi) {
    int i, j, idx, k;
    int i_s, i_f, j_s, j_f;

    double conv, pres, diff;

    double **U = sim->U;
    double **V = sim->V;

    double coef_1 = (sim->t = 1) ? -1. : -1.5;
    double coef_2 = (sim->t = 1) ? +0. : +0.5;
    double alpha = (0.5 * sim->dt) / (RE * sim->h * sim->h);

    // SOLVE FOR U

    // Implicit in X  -->  swap indices in RHS "Q" of the system
    set_system_ux(sim, adi);
    k = sim->nx + 1;

    for (int block = 0; block < 4; block++) {
        i_s = sim->i_start[block];
        i_f = sim->i_final[block];
        j_s = sim->j_start[block];
        j_f = sim->j_final[block];

        for (i = i_s; i < i_f; i++) {
            for (j = j_s; j < j_f; j++) {
                conv = coef_1 * sim->HX[i][j] + coef_2 * sim->HX_[i][j]; // swap indices since system tridiagonal with X direction
                pres = -(sim->P[i][j-1] - sim->P[i-1][j-1]) / sim->h;
                diff = U[i][j-1] - 2.*U[i][j] + U[i][j+1];

                idx = j*k + i;
                adi->a[idx] = -alpha; adi->b[idx] = 1. + 2.*alpha; adi->c[idx] = -alpha;
                adi->Q_ux[j][i] = sim->U[i][j] + 0.5*sim->dt * (conv + pres) + alpha * diff;
            }
        }
    }
    solve_thomas(sim->size_u, adi->a, adi->b, adi->c, adi->q);

    // Implicit in Y (SEEMS TO FAIL !!!)
    set_system_uy(sim, adi);
    k = sim->ny + 2;

    for (int block = 0; block < 4; block++) {
        i_s = sim->i_start[block];
        i_f = sim->i_final[block];
        j_s = sim->j_start[block];
        j_f = sim->j_final[block];

        for (i = i_s; i < i_f; i++) {
            for (j = j_s; j < j_f; j++) {
                conv = coef_1 * sim->HX[i][j] + coef_2 * sim->HX_[i][j]; // swap indices since system tridiagonal with X direction
                pres = -(sim->P[i][j-1] - sim->P[i-1][j-1]) / sim->h;
                diff = adi->Q_ux[j][i-1] - 2.*adi->Q_ux[j][i] + adi->Q_ux[j][i+1];

                idx = i*k + j;
                adi->a[idx] = -alpha; adi->b[idx] = 1. + 2.*alpha; adi->c[idx] = -alpha;
                sim->US[i][j] = adi->Q_ux[j][i] + 0.5*sim->dt * (conv + pres) + alpha * diff;
            }
        }
    }
    solve_thomas(sim->size_u, adi->a, adi->b, adi->c, sim->US[0]);


    // SOLVE FOR V

    // Implicit in X  -->  swap indices in RHS "q" of the system
    set_system_vx(sim, adi);
    k = sim->nx + 2;

    for (int block = 4; block < 8; block++) {
        i_s = sim->i_start[block];
        i_f = sim->i_final[block];
        j_s = sim->j_start[block];
        j_f = sim->j_final[block];

        for (i = i_s; i < i_f; i++) {
            for (j = j_s; j < j_f; j++) {
                conv = coef_1 * sim->HY[i][j] + coef_2 * sim->HY_[i][j]; // swap indices since system tridiagonal with X direction
                pres = -(sim->P[i-1][j] - sim->P[i-1][j-1]) / sim->h;
                diff = V[i][j-1] - 2.*V[i][j] + V[i][j+1];

                idx = j*k + i;
                adi->a[idx] = -alpha; adi->b[idx] = 1. + 2.*alpha; adi->c[idx] = -alpha;
                adi->Q_vx[j][i] = sim->V[i][j] + 0.5*sim->dt * (conv + pres) + alpha * diff;
            }
        }
    }
    solve_thomas(sim->size_v, adi->a, adi->b, adi->c, adi->q);

    // Implicit in Y
    set_system_vy(sim, adi);
    k = sim->ny + 1;

    for (int block = 4; block < 8; block++) {
        i_s = sim->i_start[block];
        i_f = sim->i_final[block];
        j_s = sim->j_start[block];
        j_f = sim->j_final[block];

        for (i = i_s; i < i_f; i++) {
            for (j = j_s; j < j_f; j++) {
                conv = coef_1 * sim->HY[i][j] + coef_2 * sim->HY_[i][j]; // swap indices since system tridiagonal with X direction
                pres = -(sim->P[i-1][j] - sim->P[i-1][j-1]) / sim->h;
                diff = adi->Q_vx[j][i-1] - 2.*adi->Q_vx[j][i] + adi->Q_vx[j][i+1];
                
                idx = i*k + j;
                adi->a[idx] = -alpha; adi->b[idx] = 1. + 2.*alpha; adi->c[idx] = -alpha;
                sim->VS[i][j] = adi->q[idx] + 0.5*sim->dt * (conv + pres) + alpha * diff;
            }
        }
    }
    solve_thomas(sim->size_v, adi->a, adi->b, adi->c, sim->VS[0]);

    // check_adi(sim);
}


void free_adi_solver(ADI_data *adi) {
    free(adi->a);
    free(adi);
}
