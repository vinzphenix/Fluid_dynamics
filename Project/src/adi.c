#include "adi.h"
#define MAX(x, y) (((x) > (y)) ? (x) : (y))

void init_adi_solver(Sim_data *sim, ADI_data *adi) {
#if USE_ADI
    int j, size, size2;
    size = sim->size_T;  // bigger than size_u and size_v
    adi->size = size;

    adi->a = (double *)malloc(4 * size * sizeof(double));
    adi->b = adi->a + size;
    adi->c = adi->b + size;
    adi->q = adi->c + size;

    size2 = (sim->ny + 2) + (sim->ny + 1) + (sim->ny + 2);

    // they point all 3 to the same array q, but at different places
    adi->Q_ux = (double **)malloc(size2 * sizeof(double *));
    adi->Q_vx = adi->Q_ux + (sim->ny + 2);
    adi->Q_Tx = adi->Q_vx + (sim->ny + 1);

    for (j = 0; j < sim->ny + 2; j++) adi->Q_ux[j] = adi->q + j * (sim->nx + 1);
    for (j = 0; j < sim->ny + 1; j++) adi->Q_vx[j] = adi->q + j * (sim->nx + 2);
    for (j = 0; j < sim->ny + 2; j++) adi->Q_Tx[j] = adi->q + j * (sim->nx + 2);
#endif
}

void test_thomas(int size, double *a, double *b, double *c, double *sol, double *rhs) {
    int i;
    double res, err;
    double tol = 1e-4;

    i = 0;
    res = b[i] * sol[i] + c[i] * sol[i + 1];
    err = fabs(res - rhs[i]);
    if (err > tol)
        printf("Error at i = %5d   ->   ERROR = %le \t rhs[i] = %lf  >< rhs[i] = %lf\n", 0, err, res, rhs[0]);

    for (i = 1; i < size - 1; i++) {
        res = a[i] * sol[i - 1] + b[i] * sol[i] + c[i] * sol[i + 1];
        err = fabs(res - rhs[i]);
        if (err > tol)
            printf("Error at i = %5d   ->   ERROR = %le \t rhs[i] = %lf  >< rhs[i] = %lf\n", i, err, res, rhs[i]);
    }

    res = a[i] * sol[i - 1] + b[i] * sol[i];
    err = fabs(res - rhs[i]);
    if (err > tol)
        printf("Error at i = %5d   ->   ERROR = %le \t rhs[i] = %le  >< rhs[i] = %le\n", i, err, res, rhs[i]);
}

#define DEBUG 0
int count = 0;
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
        w = a[i] / b[i - 1];
        b[i] -= w * c[i - 1];
        sol[i] -= w * sol[i - 1];
    }

    sol[size - 1] /= b[size - 1];
    for (i = size - 2; i >= 0; i--) {
        sol[i] = (sol[i] - c[i] * sol[i + 1]) / b[i];
    }
}

void set_system_ux(Sim_data *sim, ADI_data *adi, int k) {
    /*
     These 4 functions set_system_xx initialize the tridiagonal matrix and their corresponding RHS.
     All the nodes on the boundaries, or inside the rectangle are set to u^n.
    */

    int i, j, idx;

    // set identity for upper and lower boundaries
    for (i = 0; i < sim->nx + 1; i++) {
        for (j = 0; j < sim->ny + 2; j += sim->ny + 1) {
            idx = j * k + i;
            adi->a[idx] = 0.;
            adi->b[idx] = 1.;
            adi->c[idx] = 0.;
            adi->Q_ux[j][i] = sim->U[i][j];
        }
    }

    // set identity for left and right boundaries
    for (j = 1; j < sim->ny + 1; j++) {
        for (i = 0; i < sim->nx + 1; i += sim->nx) {
            idx = j * k + i;
            adi->a[idx] = 0.;
            adi->b[idx] = 1.;
            adi->c[idx] = 0.;
            adi->Q_ux[j][i] = sim->U[i][j];
        }
    }

    // set identity for rectangle
    for (i = sim->n * D_IN; i < (D_IN + LBOX) * sim->n + 1; i++) {
        for (j = D_BOT * sim->n + 1; j < (D_BOT + 1) * sim->n + 1; j++) {
            idx = j * k + i;
            adi->a[idx] = 0.;
            adi->b[idx] = 1.;
            adi->c[idx] = 0.;
            adi->Q_ux[j][i] = sim->U[i][j];
        }
    }
}

void set_system_uy(Sim_data *sim, ADI_data *adi, int k) {
    int i, j, idx;

    // set identity for upper and lower boundaries
    for (i = 0; i < sim->nx + 1; i++) {
        for (j = 0; j < sim->ny + 2; j += sim->ny + 1) {
            idx = i * k + j;
            adi->a[idx] = 0.;
            adi->b[idx] = 1.;
            adi->c[idx] = 0.;
            sim->US[i][j] = sim->U[i][j];
        }
    }

    // set identity for left and right boundaries
    for (j = 1; j < sim->ny + 1; j++) {
        for (i = 0; i < sim->nx + 1; i += sim->nx) {
            idx = i * k + j;
            adi->a[idx] = 0.;
            adi->b[idx] = 1.;
            adi->c[idx] = 0.;
            sim->US[i][j] = sim->U[i][j];
        }
    }

    // set identity for rectangle
    for (i = sim->n * D_IN; i < (D_IN + LBOX) * sim->n + 1; i++) {
        for (j = D_BOT * sim->n + 1; j < (D_BOT + 1) * sim->n + 1; j++) {
            idx = i * k + j;
            adi->a[idx] = 0.;
            adi->b[idx] = 1.;
            adi->c[idx] = 0.;
            sim->US[i][j] = sim->U[i][j];  // sim->uMesh;
        }
    }
}

void set_system_vx(Sim_data *sim, ADI_data *adi, int k) {
    int i, j, idx;

    // set identity for upper and lower boundaries
    for (i = 0; i < sim->nx + 2; i++) {
        for (j = 0; j < sim->ny + 1; j += sim->ny) {
            idx = j * k + i;
            adi->a[idx] = 0.;
            adi->b[idx] = 1.;
            adi->c[idx] = 0.;
            adi->Q_vx[j][i] = sim->V[i][j];
        }
    }

    // set identity for left and right boundaries
    for (j = 1; j < sim->ny; j++) {
        for (i = 0; i < sim->nx + 2; i += sim->nx + 1) {
            idx = j * k + i;
            adi->a[idx] = 0.;
            adi->b[idx] = 1.;
            adi->c[idx] = 0.;
            adi->Q_vx[j][i] = sim->V[i][j];
        }
    }

    // set identity for rectangle
    for (i = sim->n * D_IN + 1; i < (D_IN + LBOX) * sim->n + 1; i++) {
        for (j = D_BOT * sim->n; j < (D_BOT + 1) * sim->n + 1; j++) {
            idx = j * k + i;
            adi->a[idx] = 0.;
            adi->b[idx] = 1.;
            adi->c[idx] = 0.;
            adi->Q_vx[j][i] = sim->V[i][j];  // sim->vMesh;
        }
    }
}

void set_system_vy(Sim_data *sim, ADI_data *adi, int k) {
    int i, j, idx;

    // set identity for upper and lower boundaries
    for (i = 0; i < sim->nx + 2; i++) {
        for (j = 0; j < sim->ny + 1; j += sim->ny) {
            idx = i * k + j;
            adi->a[idx] = 0.;
            adi->b[idx] = 1.;
            adi->c[idx] = 0.;
            sim->VS[i][j] = sim->V[i][j];
        }
    }

    // set identity for left and right boundaries
    for (j = 1; j < sim->ny; j++) {
        for (i = 0; i < sim->nx + 2; i += sim->nx + 1) {
            idx = i * k + j;
            adi->a[idx] = 0.;
            adi->b[idx] = 1.;
            adi->c[idx] = 0.;
            sim->VS[i][j] = sim->V[i][j];
        }
    }

    // set identity for rectangle
    for (i = sim->n * D_IN + 1; i < (D_IN + LBOX) * sim->n + 1; i++) {
        for (j = D_BOT * sim->n; j < (D_BOT + 1) * sim->n + 1; j++) {
            idx = i * k + j;
            adi->a[idx] = 0.;
            adi->b[idx] = 1.;
            adi->c[idx] = 0.;
            sim->VS[i][j] = sim->V[i][j];  // sim->vMesh;
        }
    }
}

void predictor_step_u_adi(Sim_data *sim, ADI_data *adi) {
    /*
     This function updates the field u, in the first step of the LEE & KIM algorithm:
     (u* - u^n) / dt = -0.5 (3 H^n - H^(n-1)) - dP^n/dx + nu/2 (Lapl u* + Lapl u^n)
     Since the scheme is implicit in u*, we use an ADI solver to get u*.
    */

    int i, j, idx, k;
    int block, i_s, i_f, j_s, j_f;

    double **U = sim->U;

    double coef_1 = (sim->first_iteration) ? -1. : -1.5;
    double coef_2 = (sim->first_iteration) ? +0. : +0.5;
    double r_half = 0.5 * (sim->dt) / (RE * sim->h * sim->h);
    double conv, pres, diff;

    // Solve system implicit in X (reversed indices for Q_ux)
    k = sim->nx + 1;
    set_mesh_velocity(sim, sim->tnow);
    set_system_ux(sim, adi, k);

    for (block = 0; block < 4; block++) {
        i_s = sim->i_start[block];
        i_f = sim->i_final[block];
        j_s = sim->j_start[block];
        j_f = sim->j_final[block];

        for (i = i_s; i < i_f; i++) {
            for (j = j_s; j < j_f; j++) {
                conv = coef_1 * sim->HX[i][j] + coef_2 * sim->HX_[i][j];     // Convection
                pres = -(sim->P[i][j - 1] - sim->P[i - 1][j - 1]) * sim->n;  // Pressure
                diff = r_half * (U[i][j - 1] - 2. * U[i][j] + U[i][j + 1]);  // Y diffusion (explicit)

                idx = j * k + i;  // swap indices since system tridiagonal with X direction
                adi->a[idx] = -r_half;
                adi->b[idx] = 1. + 2. * r_half;
                adi->c[idx] = -r_half;
                adi->Q_ux[j][i] = sim->U[i][j] + diff + 0.5 * sim->dt * (conv + pres);
            }
        }
    }
    solve_thomas(sim->size_u, adi->a, adi->b, adi->c, adi->q);

    // Solve system implicit in Y
    k = sim->ny + 2;
    set_mesh_velocity(sim, sim->tnow + sim->dt);
    set_system_uy(sim, adi, k);

    for (block = 0; block < 4; block++) {
        i_s = sim->i_start[block];
        i_f = sim->i_final[block];
        j_s = sim->j_start[block];
        j_f = sim->j_final[block];

        for (i = i_s; i < i_f; i++) {
            for (j = j_s; j < j_f; j++) {
                // pretty dumb to compute it again, but otherwise an other array would be needed
                conv = coef_1 * sim->HX[i][j] + coef_2 * sim->HX_[i][j];
                pres = -(sim->P[i][j - 1] - sim->P[i - 1][j - 1]) * sim->n;  // equivalent to divide by dx
                diff = r_half * (adi->Q_ux[j][i - 1] - 2. * adi->Q_ux[j][i] + adi->Q_ux[j][i + 1]);

                idx = i * k + j;
                adi->a[idx] = -r_half;
                adi->b[idx] = 1. + 2. * r_half;
                adi->c[idx] = -r_half;
                sim->US[i][j] = adi->Q_ux[j][i] + diff + 0.5 * sim->dt * (conv + pres);  // US is the RHS
            }
        }
    }
    solve_thomas(sim->size_u, adi->a, adi->b, adi->c, sim->US[0]);
}

void predictor_step_v_adi(Sim_data *sim, ADI_data *adi) {
    /*
     This function updates the field v, in the first step of the LEE & KIM algorithm:
     (v* - v^n) / dt = -0.5 (3 H^n - H^(n-1)) - dP^n/dy + nu/2 (Lapl v* + Lapl v^n)
     Since the scheme is implicit in u*, we use an ADI solver to get u*.
    */

    int i, j, idx, k;
    int block, i_s, i_f, j_s, j_f;

    double **V = sim->V;

    double coef_1 = (sim->first_iteration) ? -1. : -1.5;
    double coef_2 = (sim->first_iteration) ? +0. : +0.5;
    double r_half = 0.5 * sim->dt / (RE * sim->h * sim->h);
    double conv, pres, diff;
#if TEMP_MODE
    double beta = 0.5 * GR / (RE * RE);
#endif

    // Solve system implicit in X (reversed indices for Q_ux)
    k = sim->nx + 2;
    set_mesh_velocity(sim, sim->tnow);
    set_system_vx(sim, adi, k);

    for (block = 4; block < 8; block++) {
        i_s = sim->i_start[block];
        i_f = sim->i_final[block];
        j_s = sim->j_start[block];
        j_f = sim->j_final[block];

        for (i = i_s; i < i_f; i++) {
            for (j = j_s; j < j_f; j++) {
                // swap indices since system tridiagonal with X direction
                conv = coef_1 * sim->HY[i][j] + coef_2 * sim->HY_[i][j];
                pres = -(sim->P[i - 1][j] - sim->P[i - 1][j - 1]) * sim->n;
                diff = r_half * (V[i][j - 1] - 2. * V[i][j] + V[i][j + 1]);
#if TEMP_MODE
                conv += beta * (sim->T[i][j + 1] + sim->T[i][j]);  // add the boyancy term to the RHS
#endif

                idx = j * k + i;
                adi->a[idx] = -r_half;
                adi->b[idx] = 1. + 2. * r_half;
                adi->c[idx] = -r_half;
                adi->Q_vx[j][i] = sim->V[i][j] + diff + 0.5 * sim->dt * (conv + pres);
            }
        }
    }
    solve_thomas(sim->size_v, adi->a, adi->b, adi->c, adi->q);

    // Solve system implicit in Y
    k = sim->ny + 1;
    set_mesh_velocity(sim, sim->tnow + sim->dt);
    set_system_vy(sim, adi, k);

    for (block = 4; block < 8; block++) {
        i_s = sim->i_start[block];
        i_f = sim->i_final[block];
        j_s = sim->j_start[block];
        j_f = sim->j_final[block];

        for (i = i_s; i < i_f; i++) {
            for (j = j_s; j < j_f; j++) {
                // swap indices since system tridiagonal with X direction
                conv = coef_1 * sim->HY[i][j] + coef_2 * sim->HY_[i][j];
                pres = -(sim->P[i - 1][j] - sim->P[i - 1][j - 1]) * sim->n;
                diff = r_half * (adi->Q_vx[j][i - 1] - 2. * adi->Q_vx[j][i] + adi->Q_vx[j][i + 1]);
#if TEMP_MODE
                conv += beta * (sim->T[i][j + 1] + sim->T[i][j]);  // add the boyancy term to the RHS
#endif

                idx = i * k + j;
                adi->a[idx] = -r_half;
                adi->b[idx] = 1. + 2. * r_half;
                adi->c[idx] = -r_half;
                sim->VS[i][j] = adi->Q_vx[j][i] + diff + 0.5 * sim->dt * (conv + pres);
            }
        }
    }
    solve_thomas(sim->size_v, adi->a, adi->b, adi->c, sim->VS[0]);
}

void free_adi_solver(ADI_data *adi) {
#if USE_ADI
    free(adi->a);
    free(adi->Q_ux);
#endif
    free(adi);
}
