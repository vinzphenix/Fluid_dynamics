#include "adi.h"
#define MAX(x, y) (((x) > (y)) ? (x) : (y))


void init_adi_solver(Sim_data *sim, ADI_data *adi) {
    int size = MAX(sim->size_u, sim->size_v);
    
    adi->a = (double *)calloc(4 * size, sizeof(double));
    adi->b = adi->a + size;
    adi->c = adi->b + size;
    adi->q = adi->c + size;

    /*adi->A_ux = (double **)malloc((8 * (sim->nx + 1) + 8 * (sim->nx + 2)) * sizeof(double *));
    adi->B_ux = adi->A_ux + 1 * (sim->nx + 1);
    adi->C_ux = adi->A_ux + 2 * (sim->nx + 1);
    adi->Q_ux = adi->A_ux + 3 * (sim->nx + 1);

    for (int i = 0; i < sim->nx + 1; i++) {
        adi->A_ux[i] = adi->a + i * (sim->ny + 2);
        adi->B_ux[i] = adi->b + i * (sim->ny + 2);
        adi->C_ux[i] = adi->c + i * (sim->ny + 2);
        adi->Q_ux[i] = adi->q + i * (sim->ny + 2);
    }*/
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

    // Implicit in X  -->  swap indices in RHS "q" of the system
    reset_tridiagonal(adi, sim->size_u);
    k = sim->nx + 1;
    // set us = u on boundaries
    /*for (j = 1; j < sim->ny+1; j++) {
        adi->q[j*k + 0        ] = sim->U[0      ][j];
        adi->q[j*k + (sim->nx)] = sim->U[sim->nx][j];
    }
    for (j = sim->j_final[1]; j < sim->j_start[2]; j++) {
        adi->q[j*k + sim->i_start[1]  ] = sim->U[sim->i_start[1]  ][j];
        adi->q[j*k + sim->i_final[1]-1] = sim->U[sim->i_final[1]-1][j];
    }*/
    for (i = 0; i < sim->nx + 1; i++)
        for (j = 0; j < sim->ny+2; j++)
            adi->q[j*k + i] = sim->U[i][j];

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
                adi->q[idx] = sim->U[i][j] + 0.5*sim->dt * (conv + pres) + alpha * diff;
                adi->a[idx] = -alpha;
                adi->b[idx] = 1. + 2.*alpha;
                adi->c[idx] = -alpha;
            }
        }
    }
    // for (i = (D_BOT)*sim->n * (sim->nx+1); i < ((D_BOT+1)*sim->n+2) * (sim->nx+1); i++) 
    //     PRINTF("i = %d\t a = %lf\t b = %lf\t c = %lf\t q = %lf\n", i, adi->a[i], adi->b[i], adi->c[i], adi->q[i]);
    solve_thomas(sim->size_u, adi->a, adi->b, adi->c, adi->q);
    // for (i = 0; i < sim->nx+1; i++) {
    //     for (j = (D_BOT)*sim->n; j < ((D_BOT+1)*sim->n+2); j++) {
    //         PRINTF("i = %4d\t j = %4d\t q = %lf\t u = %lf\n", i, j, adi->q[j*k+i], sim->U[i][j]);
    //     }
    //     PRINTF("\n");
    // }

    // Implicit in Y (SEEMS TO FAIL !!!)
    reset_tridiagonal(adi, sim->size_u);
    k = sim->ny + 2;
    // set us = u on boundaries
    /*for (j = 1; j < sim->ny+1; j++) {
        sim->US[0      ][j] = sim->U[0      ][j];
        sim->US[sim->nx][j] = sim->U[sim->nx][j];
    }
    for (j = sim->j_final[1]; j < sim->j_start[2]; j++) {
        sim->US[sim->i_start[1]  ][j] = sim->U[sim->i_start[1]  ][j];
        sim->US[sim->i_final[1]-1][j] = sim->U[sim->i_final[1]-1][j];
    }*/
    for (i = 0; i < sim->nx + 1; i++)
        for (j = 0; j < sim->ny+2; j++) 
            sim->US[i][j] = sim->U[i][j];

    for (int block = 0; block < 4; block++) {
        i_s = sim->i_start[block];
        i_f = sim->i_final[block];
        j_s = sim->j_start[block];
        j_f = sim->j_final[block];

        for (i = i_s; i < i_f; i++) {
            for (j = j_s; j < j_f; j++) {
                idx = j*(sim->nx + 1) + i;  // index for q
                conv = coef_1 * sim->HX[i][j] + coef_2 * sim->HX_[i][j]; // swap indices since system tridiagonal with X direction
                pres = -(sim->P[i][j-1] - sim->P[i-1][j-1]) / sim->h;
                diff = adi->q[idx-1] - 2.*adi->q[idx] + adi->q[idx+1];

                sim->US[i][j] = adi->q[idx] + 0.5*sim->dt * (conv + pres) + alpha * diff;
                idx = i*k + j;
                adi->a[idx] = -alpha;
                adi->b[idx] = 1. + 2.*alpha;
                adi->c[idx] = -alpha;
            }
        }
    }
    for (i = 0; i < sim->nx+1; i++) {
        for (j = (D_BOT)*sim->n; j < ((D_BOT+1)*sim->n + 2); j++) {
            PRINTF("i = %4d\t j = %4d\t u = %lf\t q = %lf\t us = %lf\n", i, j, sim->U[i][j], adi->q[j*(sim->nx+1)+i], sim->US[i][j]);
        }
        PRINTF("\n");
    }
    solve_thomas(sim->size_u, adi->a, adi->b, adi->c, sim->US[0]);
    // for (j = (D_BOT)*sim->n; j < ((D_BOT+1)*sim->n+2); j++)
    //     for (i = 0; i < sim->nx+1; i++)
    //         PRINTF("i = %4d\t j = %4d\t q = %lf\t u = %lf\n", i, j, sim->US[i][j], sim->U[i][j]);


    // SOLVE FOR V

    // Implicit in X  -->  swap indices in RHS "q" of the system
    reset_tridiagonal(adi, sim->size_v);
    k = sim->nx + 2;
    // set vs = v on boundaries
    /*for (i = 1; i < sim->nx+1; i++) {
        adi->q[0        *k + i] = sim->V[i][0      ];
        adi->q[(sim->ny)*k + i] = sim->V[i][sim->ny];
    }
    for (i = sim->i_start[5]; i < sim->i_start[7]; i++) {
        adi->q[ sim->j_final[5]   *k + i] = sim->V[i][sim->j_final[5]  ];
        adi->q[(sim->j_start[6]-1)*k + i] = sim->V[i][sim->j_start[6]-1];
    }*/
    for (i = 0; i < sim->nx + 2; i++)
        for (j = 0; j < sim->ny+1; j++)
            adi->q[j*k + i] = sim->V[i][j];

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
                adi->q[idx] = sim->V[i][j] + 0.5*sim->dt * (conv + pres) + alpha * diff;
                adi->a[idx] = -alpha;
                adi->b[idx] = 1. + 2.*alpha;
                adi->c[idx] = -alpha;
            }
        }
    }
    solve_thomas(sim->size_v, adi->a, adi->b, adi->c, adi->q);

    // Implicit in Y
    reset_tridiagonal(adi, sim->size_v);
    k = sim->ny + 1;
    /*for (i = 1; i < sim->nx+1; i++) {
        sim->VS[i][0      ] = sim->V[i][0      ];
        sim->VS[i][sim->ny] = sim->V[i][sim->ny];
    }
    for (i = sim->i_start[5]; i < sim->i_start[7]; i++) {
        sim->VS[i][sim->j_final[5]  ] = sim->V[i][sim->j_final[5]  ];
        sim->VS[i][sim->j_start[6]-1] = sim->V[i][sim->j_start[6]-1];
    }*/
    for (i = 0; i < sim->nx + 2; i++)
        for (j = 0; j < sim->ny+1; j++) 
            sim->VS[i][j] = sim->V[i][j];

    for (int block = 4; block < 8; block++) {
        i_s = sim->i_start[block];
        i_f = sim->i_final[block];
        j_s = sim->j_start[block];
        j_f = sim->j_final[block];

        for (i = i_s; i < i_f; i++) {
            for (j = j_s; j < j_f; j++) {
                idx = j*(sim->nx + 2) + i;  // index for q
                conv = coef_1 * sim->HY[i][j] + coef_2 * sim->HY_[i][j]; // swap indices since system tridiagonal with X direction
                pres = -(sim->P[i-1][j] - sim->P[i-1][j-1]) / sim->h;
                diff = adi->q[idx-1] - 2.*adi->q[idx] + adi->q[idx+1];
                
                sim->VS[i][j] = adi->q[idx] + 0.5*sim->dt * (conv + pres) + alpha * diff;
                idx = i*k + j;
                adi->a[idx] = -alpha;
                adi->b[idx] = 1. + 2.*alpha;
                adi->c[idx] = -alpha;
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
