#include "poisson.h"
#include "project.h"
#include "adi.h"

char *myPath = "./data/";
char filename_params[50];
char filename_u[50];
char filename_v[50];
char filename_p[50];
char filename_w[50];


void init_Sim_data(Sim_data *sim) {
    sprintf(filename_params, "%ssimu_params.txt", myPath);
    sprintf(filename_u, "%ssimu_u.txt", myPath);
    sprintf(filename_v, "%ssimu_v.txt", myPath);
    sprintf(filename_p, "%ssimu_p.txt", myPath);
    sprintf(filename_w, "%ssimu_w.txt", myPath);


    // sim->h = CFL / (FOURIER * RE * 5);
    // sim->dt = FOURIER * RE * (sim->h) * (sim->h);

    // sim->nt = ceil(TEND / sim->dt);
    // sim->nx = (double) L_ / sim->h;
    // sim->ny = (double) H_ / sim->h;

#if TEST_POISSON
    sim->n = 4;
#else
    sim->n = N;  // 1. / h
#endif

    sim->nt = NT;
    sim->nx = L_ * sim->n;
    sim->ny = H_ * sim->n;

    // sim->dt = TEND / (double) sim->nt;
    sim->dt = 0.001;
    sim->h = 1. / ((double) sim->n);
    
    sim->uMesh = 0.;
    sim->vMesh = 0.;

    int size_u = (sim->nx + 1) * (sim->ny + 2);
    int size_v = (sim->nx + 2) * (sim->ny + 1);
    int size_p = (sim->nx + 0) * (sim->ny + 0);
    sim->size_u = size_u;
    sim->size_v = size_v;
    sim->size_p = size_p;
    

    int n_vectors = 4;

    // Allocate u and set matrix access
    sim->u_data = (double *)calloc(n_vectors * size_u, sizeof(double));

    sim->U = (double **)malloc(4*(sim->nx + 1) * sizeof(double *));
    sim->US  = sim->U + 1 * (sim->nx + 1);
    sim->HX  = sim->U + 2 * (sim->nx + 1);
    sim->HX_ = sim->U + 3 * (sim->nx + 1);

    for (int i = 0; i < sim->nx + 1; i++) {
        sim->U[i]   = sim->u_data + 0 * size_u + i * (sim->ny + 2);
        sim->US[i]  = sim->u_data + 1 * size_u + i * (sim->ny + 2);
        sim->HX[i]  = sim->u_data + 2 * size_u + i * (sim->ny + 2);
        sim->HX_[i] = sim->u_data + 3 * size_u + i * (sim->ny + 2);
    }


    // Allocate v and set matrix access
    sim->v_data = (double *)calloc(n_vectors * size_v, sizeof(double));
    
    sim->V = (double **)malloc(4*(sim->nx + 2) * sizeof(double *));
    sim->VS  = sim->V + 1 * (sim->nx + 2);
    sim->HY  = sim->V + 2 * (sim->nx + 2);
    sim->HY_ = sim->V + 3 * (sim->nx + 2);

    for (int i = 0; i < sim->nx + 2; i++) {
        sim->V[i]   = sim->v_data + 0 * size_v + i * (sim->ny + 1);
        sim->VS[i]  = sim->v_data + 1 * size_v + i * (sim->ny + 1);
        sim->HY[i]  = sim->v_data + 2 * size_v + i * (sim->ny + 1);
        sim->HY_[i] = sim->v_data + 3 * size_v + i * (sim->ny + 1);
    }


    // Allocate p and set matrix access
    sim->p_data = (double *)calloc(2 * size_p, sizeof(double));

    sim->P = (double **)malloc(2 * (sim->nx) * sizeof(double *));
    sim->PHI = sim->P + (sim->nx);

    for (int i = 0; i < sim->nx; i++) {
        sim->P[i]   = sim->p_data + 0 * size_p + i * (sim->ny);
        sim->PHI[i] = sim->p_data + 1 * size_p + i * (sim->ny);
    }


    // set bounds for sub-domains
    sim->i_start[0] = 1;                          sim->i_final[0] = (D_IN * sim->n);
    sim->i_start[1] = (D_IN * sim->n);            sim->i_final[1] = (D_IN + LBOX) * sim->n + 1;
    sim->i_start[2] = (D_IN * sim->n);            sim->i_final[2] = (D_IN + LBOX) * sim->n + 1;
    sim->i_start[3] = (D_IN + LBOX) * sim->n + 1; sim->i_final[3] = sim->nx;

    sim->j_start[0] = 1;                          sim->j_final[0] = sim->ny + 1;
    sim->j_start[1] = 1;                          sim->j_final[1] = (D_BOT) * sim->n + 1;
    sim->j_start[2] = (D_BOT + 1) * sim->n + 1;   sim->j_final[2] = sim->ny + 1;
    sim->j_start[3] = 1;                          sim->j_final[3] = sim->ny + 1;

    sim->i_start[4] = sim->i_start[0];            sim->i_final[4] = sim->i_final[0] + 1;
    sim->i_start[5] = sim->i_start[1] + 1;        sim->i_final[5] = sim->i_final[1];
    sim->i_start[6] = sim->i_start[2] + 1;        sim->i_final[6] = sim->i_final[2];
    sim->i_start[7] = sim->i_start[3];            sim->i_final[7] = sim->i_final[3] + 1;

    sim->j_start[4] = sim->j_start[0];            sim->j_final[4] = sim->j_final[0] - 1;
    sim->j_start[5] = sim->j_start[1];            sim->j_final[5] = sim->j_final[1] - 1;
    sim->j_start[6] = sim->j_start[2];            sim->j_final[6] = sim->j_final[2] - 1;
    sim->j_start[7] = sim->j_start[3];            sim->j_final[7] = sim->j_final[3] - 1;

}



void init_fields(Sim_data *sim) {
    int i, j;

    int i_w_left =  D_IN * sim->n;           // on boundary
    int i_w_right = (D_IN + LBOX) * sim->n;  // on boundary
    int j_w_below = D_BOT * sim->n + 1;      // in rectangle
    int j_w_above = (D_BOT + 1) * sim->n;    // in rectangle


    // set u to Uinf at beginning (0 is good for v)
    for (i = 0; i < sim->nx+1; i++) {
        for (j = 0; j < sim->ny+2; j++) {
            if ((sim->i_final[0] <= i) && (i < sim->i_start[3]) && (sim->j_final[1] <= j) && (j < sim->j_start[2])) {
                // printf("x = %.3lf  y = %.3lf\n", i*sim->h, (j-0.5)*sim->h);
                continue;
            } else {
                sim->U[i][j] = 1.;
                sim->US[i][j] = 1.;
            }
        }
    }

    /*// set u to 0 inside rectangle (to be sure it stays constant)
    for (i = i_w_left; i < i_w_right + 1; i++) {
        for (j = j_w_below; j < j_w_above + 1; j++) {
            sim->u[i*k + j] = 0.;
        }
    }

    // set u to 0 on left and right walls
    for (j = j_w_below; j < j_w_above + 1; j++) {
        i = i_w_left;
        sim->u[i*k + j] = 0.;
        
        i = i_w_right;
        sim->u[i*k + j] = 0.;
    }*/

    // set ghost for u on upper and lower walls of rectangle
    for (i = i_w_left + 1; i < i_w_right; i++) {
        j = j_w_below;
        sim->U[i][j] = -0.2 * 11.; // since u(t=0) is 1.

        j = j_w_above;
        sim->U[i][j] = -0.2 * 11.; // since u(t=0) is 1.
    }
}


void save_fields(Sim_data *sim, int t) {
    FILE *ptr, *ptr_u, *ptr_v, *ptr_p, *ptr_w;

    if (t == 0) {
        ptr = fopen(filename_params, "w");
        fprintf(ptr, "%d %d %d\n", sim->nt / SAVE_MODULO, sim->nx, sim->ny);
        fprintf(ptr, "%lf %d %d %d %d %d %lf %lf %lf %lf %lf\n", TEND, L_, H_, LBOX, D_IN, D_BOT, sim->dt, sim->h, ALPHA, STROUHAL, sim->dt * SAVE_MODULO);
        fclose(ptr);
    }

    int i, j;
    // int k = sim->ny+2;
    // int l = sim->ny+1;
    double u_up, u_dw, v_lf, v_rg;
    
    if (t == 0) {
        ptr_u = fopen(filename_u, "w");
        ptr_v = fopen(filename_v, "w");
        ptr_w = fopen(filename_w, "w");
        ptr_p = fopen(filename_p, "w");
    } else {
        ptr_u = fopen(filename_u, "a");
        ptr_v = fopen(filename_v, "a");
        ptr_w = fopen(filename_w, "a");
        ptr_p = fopen(filename_p, "a");
    }

    // loop over corners
    for (i = 0; i < sim->nx + 1; i++) {
        for (j = 0; j < sim->ny + 1; j++) {
            
            u_up = sim->U[i][j + 1]; // u above
            u_dw = sim->U[i][j];     // u below
            v_lf = sim->V[i][j];     // v left
            v_rg = sim->V[i + 1][j]; // v right
            
            fprintf(ptr_u, FMT, 0.5 * (u_up + u_dw));                       // u average 
            fprintf(ptr_v, FMT, 0.5 * (v_lf + v_rg));                       // v average
            fprintf(ptr_w, FMT, ((v_rg - v_lf) - (u_up - u_dw)) / sim->h);  // vorticity
        }
    }

    // loop over centers
    for (i = 0; i < sim->size_p; i++) {
        fprintf(ptr_p, FMT, sim->p_data[i]);  // pressure
    }
    
    // fprintf(ptr_u, "\n");
    fclose(ptr_u);
    // fprintf(ptr_v, "\n");
    fclose(ptr_v);
    // fprintf(ptr_w, "\n");
    fclose(ptr_w);
    // fprintf(ptr_p, "\n");
    fclose(ptr_p);

}


void set_bd_conditions(Sim_data *sim, double **U, double **V) {
    int i, j;
    
    // int k = sim->ny + 2;
    // int l = sim->ny + 1;

    /*// Inflow boundary condition for u
    for (j = 1; j <= sim->ny; j++) {   // Inflow uniform profile u
        U[0][j] = 1.;   // could do a special profile
    }*/

    // Outflow boundary condition for u
    i = sim->nx;
    for (j = 1; j <= sim->ny; j++) {  // u is advected at velocity (1 - uMesh)
#if OUTFLOW_STUPID
        U[i][j] = 1.;
#else
        U[i][j] -= sim->dt / sim->h * (1. - sim->uMesh) * (sim->U[i][j] - sim->U[i-1][j]);
#endif
    }

    /*// Lateral walls condition for v
    for (i = 0; i < sim->nx + 2; i++) { // no-through flow
        V[i][0] = 0.;             // below
        V[i][sim->ny] = 0.;       // above
    }*/


    i = D_IN * sim->n;  // Left wall of rectangle (u)
    for (j = D_BOT * sim->n + 1; j <= (D_BOT + 1) * sim->n; j++)
        U[i][j] = sim->uMesh;
    
    i = (D_IN + LBOX) * sim->n; // Right wall of rectangle (u)
    for (j = D_BOT * sim->n + 1; j <= (D_BOT + 1) * sim->n; j++)
        U[i][j] = sim->uMesh;


    for (i = D_IN * sim->n + 1; i < (D_IN + LBOX) * sim->n; i++) {                
        
        j = D_BOT * sim->n;  // lower wall of the rectangle
        V[i][j] = sim->vMesh;

        j = (D_BOT + 1) * sim->n;  // upper wall of the rectangle
        V[i][j] = sim->vMesh;
        // printf("x = %.3lf,  y = %.3lf,  v[%d, %d] = %.3lf\n", (i-0.5)*sim->h, j*sim->h, i,j, v[i, j]); fflush(stdout);
    }
}

void set_ghost_points(Sim_data *sim) {
    int i, j;

    // External walls (v = 0 and w = 0 --> du/dy = 0)
    for (i = 0; i < sim->nx + 1; i++) {                // zero voriticty
        sim->U[i][0] = sim->U[i][1];                   // below
        sim->U[i][sim->ny + 1] = sim->U[i][sim->ny];   // above
    }

    // Inflow boundary condition for v
    for (j = 1; j < sim->ny; j++) {  // trick with ghost for zero v at inflow
        sim->V[0][j] = -0.2 * (sim->V[3][j] - 5. * sim->V[2][j] + 15. * sim->V[1][j] - 16. * 0.);  // 0 or vMesh ?
    }
    
    // Outflow boundary condition for v

#if OUTFLOW_STUPID
    for (j = 1; j < sim->ny; j++) {
        sim->V[sim->nx+1][j] = 0.;
    }

#else
    double w_last, w_left;
    i = sim->nx;
    for (j = 1; j < sim->ny; j++) {
        w_last = sim->V[i+1][j] - sim->V[i  ][j] - sim->U[i  ][j+1] + sim->U[i  ][j];
        w_left = sim->V[i  ][j] - sim->V[i-1][j] - sim->U[i-1][j+1] + sim->U[i-1][j];
        sim->V[i+1][j] = (sim->V[i][j] + sim->U[i][j+1] - sim->U[i][j]) \
                           + w_last - sim->dt / sim->h * (1 - sim->uMesh) * (w_last - w_left);
        // [ v_(i+0.5,j) - v_(i-0.5,j) - u_(i,j+0.5) + u(i,j-0.5) ] / h = w_ij + dt uc / h * [ w_ij - w_(i-1, j) ]
    }
#endif

    // Left wall of rectangle
    i = D_IN * sim->n + 1;
    for (j = D_BOT * sim->n + 1; j < (D_BOT + 1) * sim->n; j++)
        sim->V[i][j] = -0.2 * (sim->V[i-3][j] - 5.*sim->V[i-2][j] + 15.*sim->V[i-1][j] - 16.*sim->vMesh);

    // Right wall of rectangle
    i = (D_IN + LBOX) * sim->n;
    for (j = D_BOT * sim->n + 1; j < (D_BOT + 1) * sim->n; j++)
        sim->V[i][j] = -0.2 * (sim->V[i+3][j] - 5.*sim->V[i+2][j] + 15.*sim->V[i+1][j] - 16.*sim->vMesh);

    // Upper and lower walls of rectangle
    for (i = D_IN * sim->n + 1; i < (D_IN + LBOX) * sim->n; i++) {
    
        j = D_BOT * sim->n + 1;  // below
        sim->U[i][j] = -0.2 * (sim->U[i][j-3] - 5.*sim->U[i][j-2] + 15.*sim->U[i][j-1] - 16*sim->uMesh);

        j = (D_BOT + 1) * sim->n;  // above
        sim->U[i][j] = -0.2 * (sim->U[i][j+3] - 5.*sim->U[i][j+2] + 15.*sim->U[i][j+1] - 16*sim->uMesh);
    }
}


void compute_convection(Sim_data *sim) {
    int i, j;

    double **U = sim->U;
    double **V = sim->V;
    double uMesh = sim->uMesh;
    double vMesh = sim->vMesh;

    double **H = sim->HX;
    /*for (i = 1; i < sim->nx; i++) {
        for (j = 1; j <= sim->ny; j++) {*/

    int i_s, i_f, j_s, j_f;

    for (int block = 0; block < 4; block++) {
        i_s = sim->i_start[block];
        i_f = sim->i_final[block];
        j_s = sim->j_start[block];
        j_f = sim->j_final[block];

        for (i = i_s;  i < i_f; i++) {
            for (j = j_s; j < j_f; j++) {

                H[i][j] = 0.;

                // Advective form for u field
    #if         (CONVECTION_MODE == 0) || (CONVECTION_MODE == 2)
                H[i][j] += (U[i+1][j  ] + U[i  ][j  ] - 2.*uMesh) * (U[i+1][j  ] - U[i  ][j  ])\
                         + (U[i  ][j  ] + U[i-1][j  ] - 2.*uMesh) * (U[i  ][j  ] - U[i-1][j  ])\
                         + (V[i  ][j  ] + V[i+1][j  ] - 2.*vMesh) * (U[i  ][j+1] - U[i  ][j  ])\
                         + (V[i  ][j-1] + V[i+1][j-1] - 2.*vMesh) * (U[i  ][j  ] - U[i  ][j-1]);
    #endif
    #if         (CONVECTION_MODE == 1) || (CONVECTION_MODE == 2)
                // Divergence form for u field
                H[i][j] += (U[i+1][j  ] + U[i][j]) * (U[i+1][j  ] + U[i  ][j  ] - 2.*uMesh)\
                         - (U[i-1][j  ] + U[i][j]) * (U[i-1][j  ] + U[i  ][j  ] - 2.*uMesh)\
                         + (U[i  ][j+1] + U[i][j]) * (V[i  ][j  ] + V[i+1][j  ] - 2.*vMesh)\
                         - (U[i  ][j-1] + U[i][j]) * (V[i  ][j-1] + V[i+1][j-1] - 2.*vMesh);
    #endif
                H[i][j] /= (4. * sim->h);  // possible since dx = dy = h; factor "2" since avg of Adv / Div

    #if         (CONVECTION_MODE == 2)
                H[i][j] /= 2.;
    #endif
            }
        }
    }

    H = sim->HY;
    for (int block = 4; block < 8; block++) {
        i_s = sim->i_start[block];
        i_f = sim->i_final[block];
        j_s = sim->j_start[block];
        j_f = sim->j_final[block];
    
        /*for (i = 1; i <= sim->nx; i++) {
            for (j = 1; j < sim->ny; j++) {*/

        for (i = i_s;  i < i_f; i++) {
            for (j = j_s; j < j_f; j++) {
                
                H[i][j] = 0.;

    #if         (CONVECTION_MODE == 0) || (CONVECTION_MODE == 2)
                // Advective form for v field
                H[i][j] += (U[i  ][j+1] + U[i  ][j  ] - 2.*uMesh) * (V[i+1][j  ] - V[i  ][j  ])
                         + (U[i-1][j+1] + U[i-1][j  ] - 2.*uMesh) * (V[i  ][j  ] - V[i-1][j  ])
                         + (V[i  ][j+1] + V[i  ][j  ] - 2.*vMesh) * (V[i  ][j+1] - V[i  ][j  ])
                         + (V[i  ][j  ] + V[i  ][j-1] - 2.*vMesh) * (V[i  ][j  ] - V[i  ][j-1]);
    #endif

    #if         (CONVECTION_MODE == 1) || (CONVECTION_MODE == 2)
                // Divergence form for v field
                H[i][j] += (V[i+1][j  ] + V[i][j]) * (U[i  ][j+1] + U[i  ][j  ] - 2.*uMesh)
                         - (V[i-1][j  ] + V[i][j]) * (U[i-1][j+1] + U[i-1][j  ] - 2.*uMesh)
                         + (V[i  ][j+1] + V[i][j]) * (V[i  ][j+1] + V[i  ][j  ] - 2.*vMesh)
                         - (V[i  ][j-1] + V[i][j]) * (V[i  ][j-1] + V[i  ][j  ] - 2.*vMesh);
    #endif
                H[i][j] /= (4. * sim->h);

    #if         (CONVECTION_MODE == 2)
                H[i][j] /= 2.;
    #endif

            }
        }
    }
}


void check_convection(Sim_data *sim) {
    int i, j;
    double x, y;
    FILE *ptr = fopen("check_convection.txt", "w");

    for (i = 0; i < sim->nx+1; i++) {
        for (j = 0; j < sim->ny+2; j++) {
            x = i * sim->h;
            y = (j - 0.5) * sim->h;
            fprintf(ptr, "%s %le %le %le %le\n", "u", x, y, sim->U[i][j], sim->US[i][j]);
        }
    }
    for (i = 0; i < sim->nx+2; i++) {
        for (j = 0; j < sim->ny+1; j++) {
            x = (i - 0.5) * sim->h;
            y = j * sim->h;
            fprintf(ptr, "%s %le %le %le %le\n", "v", x, y, sim->V[i][j], sim->VS[i][j]);
        }
    }
    for (i = 0; i < sim->nx+1; i++) {
        for (j = 0; j < sim->ny+2; j++) {
            x = i * sim->h;
            y = (j - 0.5) * sim->h;
            fprintf(ptr, "%c %le %le %le %le\n", 'X', x, y, sim->HX[i][j], sim->HX_[i][j]);
        }
    }
    for (i = 0; i < sim->nx+2; i++) {
        for (j = 0; j < sim->ny+1; j++) {
            x = (i - 0.5) * sim->h;
            y = j * sim->h;
            fprintf(ptr, "%c %le %le %le %le\n", 'Y', x, y, sim->HY[i][j], sim->HY_[i][j]);
        }
    }
    fclose(ptr);
}

void save_debug(Sim_data *sim) {
    int i, j;
    FILE *ptr_u = fopen("fields_u.txt", "w");
    FILE *ptr_v = fopen("fields_v.txt", "w");

    for (j = sim->ny+1; j >= 0; j--) {
        for (i = 0 ; i < sim->nx+1; i++) {
            fprintf(ptr_u, "%3.0lf/%7.2lf     ", sim->U[i][j], sim->HX[i][j]);
        }
        fprintf(ptr_u, "\n");
    }

    for (j = sim->ny; j >= 0; j--) {
        for (i = 0 ; i <= sim->nx+1; i++) {
            fprintf(ptr_v, "%3.0lf/%7.2lf     ", sim->V[i][j], sim->HY[i][j]);
        }
        fprintf(ptr_v, "\n");
    }

    fclose(ptr_u);
    fclose(ptr_v);
    printf("done\n");
}


void predictor_step(Sim_data *sim) {
    int i, j;
    double conv, pres, diff;

    double **U = sim->U;
    double **V = sim->V;

    double coef_1 = (sim->t = 1) ? -1. : -1.5;
    double coef_2 = (sim->t = 1) ? +0. : +0.5;
    double alpha = 1. / (RE * sim->h * sim->h);

    // update u field
    /*for (i = 1; i < sim->nx; i++) {
        for (j = 1; j <= sim->ny; j++) {*/

    int i_s, i_f, j_s, j_f;

    for (int block = 0; block < 4; block++) {
        i_s = sim->i_start[block];
        i_f = sim->i_final[block];
        j_s = sim->j_start[block];
        j_f = sim->j_final[block];

        for (i = i_s;  i < i_f; i++) {
            for (j = j_s; j < j_f; j++) {
                            
                conv = coef_1 * sim->HX[i][j] + coef_2 * sim->HX_[i][j];  // (u du/dx + v du/dy) at t_n and t_(n-1)
                pres = -(sim->P[i][j-1] - sim->P[i-1][j-1]) / sim->h;  // dp/dx
                diff = alpha * (U[i+1][j] + U[i-1][j] - 4*U[i][j] + U[i][j+1] + U[i][j-1]);  // d2u/dx2 + d2u/dy2
                
                sim->US[i][j] = U[i][j] + sim->dt * (conv + pres + diff);
            }
        }
    }

    // update v field
    /*for (i = 1; i <= sim->nx; i++) {
        for (j = 1; j < sim->ny; j++) {*/

    for (int block = 4; block < 8; block++) {
        i_s = sim->i_start[block];
        i_f = sim->i_final[block];
        j_s = sim->j_start[block];
        j_f = sim->j_final[block];

        for (i = i_s;  i < i_f; i++) {
            for (j = j_s; j < j_f; j++) {
                
                conv = coef_1 * sim->HY[i][j] + coef_2 * sim->HY_[i][j];  // (u dv/dx + v dv/dy) at t_n and t_(n-1)
                pres = -(sim->P[i-1][j] - sim->P[i-1][j-1]) / sim->h;  // dp/dy
                diff = alpha * (V[i+1][j] + V[i-1][j] - 4*V[i][j] + V[i][j+1] + V[i][j-1]);  // d2v/dx2 + d2v/dy2

                sim->VS[i][j] = V[i][j] + sim->dt * (conv + pres + diff);
            }
        }
    }
}

void corrector_step(Sim_data *sim) {
    int i, j;
    double dphi;

    // update u field
    /*for (i = 1; i < sim->nx; i++) {
        for (j = 1; j <= sim->ny; j++) {*/

    int i_s, i_f, j_s, j_f;

    for (int block = 0; block < 4; block++) {
        i_s = sim->i_start[block];
        i_f = sim->i_final[block];
        j_s = sim->j_start[block];
        j_f = sim->j_final[block];
        
        for (i = i_s;  i < i_f; i++) {
            for (j = j_s; j < j_f; j++) {
                dphi = (sim->PHI[i][j-1] - sim->PHI[i-1][j-1]) / sim->h;
                sim->U[i][j] = sim->US[i][j] - sim->dt * dphi;  // - d phi / dx
            }
        }
    }

    // update v field
    /*for (i = 1; i <= sim->nx; i++) {
        for (j = 1; j < sim->ny; j++) {*/

    for (int block = 4; block < 8; block++) {
        i_s = sim->i_start[block];
        i_f = sim->i_final[block];
        j_s = sim->j_start[block];
        j_f = sim->j_final[block];
        
        for (i = i_s;  i < i_f; i++) {
            for (j = j_s; j < j_f; j++) {
                dphi = (sim->PHI[i-1][j] - sim->PHI[i-1][j-1]) / sim->h;
                sim->V[i][j] = sim->VS[i][j] - sim->dt * dphi; // - d phi / dy
            }
        }
    }

    // update p field
    for (i = 0, j = sim->size_p; i < sim->size_p; i++, j++) {
        sim->p_data[i] += sim->p_data[j];
    }

    /*for (i = 0; i < sim->nx; i++) {
        for (j = 0; j < sim->ny; j++) {
            sim->P[i][j] += sim->PHI[i][j];
        }
    }*/
}


void swap_next_previous(Sim_data *sim) {
    double **temp;
    
    temp = sim->HX;
    sim->HX = sim->HX_;
    sim->HX_ = temp;

    temp = sim->HY;
    sim->HY = sim->HY_;
    sim->HY_ = temp;
}


void integrate_flow(Sim_data *sim, Poisson_data *poisson, ADI_data *adi) {
    int t = 0;

    set_ghost_points(sim);
    set_bd_conditions(sim, sim->U, sim->V);
    set_bd_conditions(sim, sim->US, sim->VS);
    save_fields(sim, t);

    printf("starting ... \n"); fflush(stdout);

    for (t = 1; t <= sim->nt; t++) {
        
        sim->t = t;
        if (T_START < t * sim->dt) {
            sim->uMesh = ALPHA * 2. * M_PI * STROUHAL * sin(2. * M_PI * STROUHAL * (t * sim->dt - T_START));
            sim->vMesh = 0.;
        }
        
        // set_bd_conditions(sim, sim->U, sim->V);
        // set_ghost_points(sim);

        compute_convection(sim);

#       if USE_ADI
        predictor_step_adi(sim, adi);
#       else
        predictor_step(sim);
#       endif

        set_bd_conditions(sim, sim->US, sim->VS);

        // check_convection(sim);

        poisson_solver(sim, poisson);
        corrector_step(sim);
        set_bd_conditions(sim, sim->U, sim->V);
        set_ghost_points(sim);

        swap_next_previous(sim);

        if (t % SAVE_MODULO == 0) {
            printf("Saving results ... t = %3d\n", t);
            save_fields(sim, t);
        }
    }

    printf("\n");

}


void free_Sim_data(Sim_data *sim) {
    free(sim->u_data); free(sim->U);
    free(sim->v_data); free(sim->V);
    free(sim->p_data); free(sim->P);
    free(sim);
}


int main(int argc, char *argv[]){
    // argv : -ksp_type gmres -pc_type lu
    PetscInitialize(&argc, &argv, 0, 0);

    Sim_data *simulation = (Sim_data *)malloc(sizeof(Sim_data));
    init_Sim_data(simulation);
    init_fields(simulation);

    Poisson_data *poisson = (Poisson_data *)malloc(sizeof(Poisson_data));
    initialize_poisson_solver(simulation, poisson);

    ADI_data *adi_solver = (ADI_data *)malloc(sizeof(ADI_data));
#   if USE_ADI
    init_adi_solver(simulation, adi_solver);
#   endif


    // MAIN PROCESS
#   if TEST_POISSON    
    test_poisson(simulation, poisson);
#   else
    integrate_flow(simulation, poisson, adi_solver);
#   endif


    // Free memory
    free_Sim_data(simulation);
    free_poisson_solver(poisson);
#   if USE_ADI
    free_adi_solver(adi_solver);
#   else
    free(adi_solver);
#   endif
    PetscFinalize();
}
