#include "poisson.h"
#include "project.h"
#include "adi.h"

char *myPath = "./results";
char filename_params[50];
char filename_u[50];
char filename_v[50];
char filename_p[50];
char filename_T[50];


void init_Sim_data(Sim_data *sim, int n_input, double dt_input, double tend_input) {
    sprintf(filename_params, "%s/simu_params.txt", myPath);
    sprintf(filename_u, "%s/simu_u.txt", myPath);
    sprintf(filename_v, "%s/simu_v.txt", myPath);
    sprintf(filename_p, "%s/simu_p.txt", myPath);
    sprintf(filename_T, "%s/simu_T.txt", myPath);

    // Space discretization
    // sim->n = N_;
    sim->n = n_input;
    sim->nx = L_ * sim->n;
    sim->ny = H_ * sim->n;
    sim->h = 1. / ((double) sim->n);

    // Time discretization
    // sim->dt = DT;
    // sim->tsim = TSIM;
    sim->dt = dt_input;
    sim->tsim = tend_input;
    sim->nt = (int) ceil(sim->tsim / sim->dt);
    double dtStable = fmin(FOURIER * RE * (sim->h) * (sim->h), CFL * sim->h / U0V0);
    printf("Current \u0394t = %.5lf  vs  %.5lf = \u0394t maximal\n", sim->dt, dtStable);
    
    sim->uMesh = 0.;
    sim->vMesh = 0.;

    int size_u = (sim->nx + 1) * (sim->ny + 2);
    int size_v = (sim->nx + 2) * (sim->ny + 1);
    int size_p = (sim->nx + 0) * (sim->ny + 0);
    int size_T = (sim->nx + 2) * (sim->ny + 2);
    sim->size_u = size_u;
    sim->size_v = size_v;
    sim->size_p = size_p;
    sim->size_T = size_T;
    

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
    
    sim->V = (double **)malloc(4 * (sim->nx + 2) * sizeof(double *));
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


#   if TEMPERATURE
    // Allocate T and set matrix access
    sim->T_data = (double *)calloc(3 * size_T, sizeof(double));

    sim->T = (double **)malloc(3 * (sim->nx + 2) * sizeof(double *));
    sim->HT  = sim->T + 1 * (sim->nx + 2);
    sim->HT_ = sim->T + 2 * (sim->nx + 2);

    for (int i = 0; i < sim->nx + 2; i++) {
        sim->T[i]   = sim->T_data + 0 * size_T + i * (sim->ny + 2);
        sim->HT[i]  = sim->T_data + 1 * size_T + i * (sim->ny + 2);
        sim->HT_[i] = sim->T_data + 2 * size_T + i * (sim->ny + 2);
    }
#   endif


    // set bounds for sub-domains
    // for u
    sim->i_start[0] = 1;                          sim->i_final[0] = (D_IN * sim->n);
    sim->i_start[1] = (D_IN * sim->n);            sim->i_final[1] = (D_IN + LBOX) * sim->n + 1;
    sim->i_start[2] = (D_IN * sim->n);            sim->i_final[2] = (D_IN + LBOX) * sim->n + 1;
    sim->i_start[3] = (D_IN + LBOX) * sim->n + 1; sim->i_final[3] = sim->nx;

    sim->j_start[0] = 1;                          sim->j_final[0] = sim->ny + 1;
    sim->j_start[1] = 1;                          sim->j_final[1] = (D_BOT) * sim->n + 1;
    sim->j_start[2] = (D_BOT + 1) * sim->n + 1;   sim->j_final[2] = sim->ny + 1;
    sim->j_start[3] = 1;                          sim->j_final[3] = sim->ny + 1;

    // for v
    sim->i_start[4] = sim->i_start[0];            sim->i_final[4] = sim->i_final[0] + 1;
    sim->i_start[5] = sim->i_start[1] + 1;        sim->i_final[5] = sim->i_final[1];
    sim->i_start[6] = sim->i_start[2] + 1;        sim->i_final[6] = sim->i_final[2];
    sim->i_start[7] = sim->i_start[3];            sim->i_final[7] = sim->i_final[3] + 1;

    sim->j_start[4] = sim->j_start[0];            sim->j_final[4] = sim->j_final[0] - 1;
    sim->j_start[5] = sim->j_start[1];            sim->j_final[5] = sim->j_final[1] - 1;
    sim->j_start[6] = sim->j_start[2];            sim->j_final[6] = sim->j_final[2] - 1;
    sim->j_start[7] = sim->j_start[3];            sim->j_final[7] = sim->j_final[3] - 1;

#   if TEMPERATURE
    sim->i_start[8 ] = sim->i_start[0];            sim->i_final[8 ] = sim->i_final[0] + 1;
    sim->i_start[9 ] = sim->i_start[1] + 1;        sim->i_final[9 ] = sim->i_final[1];
    sim->i_start[10] = sim->i_start[2] + 1;        sim->i_final[10] = sim->i_final[2];     
    sim->i_start[11] = sim->i_start[3];            sim->i_final[11] = sim->i_final[3] + 1;

    sim->j_start[8 ] = sim->j_start[0];            sim->j_final[8 ] = sim->j_final[0];
    sim->j_start[9 ] = sim->j_start[1];            sim->j_final[9 ] = sim->j_final[1];
    sim->j_start[10] = sim->j_start[2];            sim->j_final[10] = sim->j_final[2];
    sim->j_start[11] = sim->j_start[3];            sim->j_final[11] = sim->j_final[3];
#   endif
}


void init_fields(Sim_data *sim) {
    int i, j;

    int i_w_left =  D_IN * sim->n;           // on boundary
    int i_w_right = (D_IN + LBOX) * sim->n;  // on boundary
    int j_w_below = D_BOT * sim->n + 1;      // in rectangle
    int j_w_above = (D_BOT + 1) * sim->n;    // in rectangle


    // set u to Uinf at beginning (0 is already good for v)
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

    // set ghost for u on upper and lower walls of rectangle
    for (i = i_w_left + 1; i < i_w_right; i++) {
        j = j_w_below;
        sim->U[i][j] = -0.2 * (sim->U[i][j-3] - 5.*sim->U[i][j-2] + 15.*sim->U[i][j-1] - 16 * 0.);

        j = j_w_above;
        sim->U[i][j] = -0.2 * (sim->U[i][j+3] - 5.*sim->U[i][j+2] + 15.*sim->U[i][j+1] - 16 * 0.);
    }

    
#   if NO_SLIP
    double eta;
    for (j = 1; j < sim->ny + 1; j++) {   // Inflow uniform profile u
        eta = (j-0.5) * sim->h / H_;
        sim->U[0][j] = 6. * (1.) * eta * (1. - eta);  // Poiseuille
        sim->US[0][j] = sim->U[0][j];                  // Poiseuille
    }
    for (i = 0; i < sim->nx + 1; i++) {
        j = 0;
        sim->U[i][j] = -0.2 * (sim->U[i][j+3] - 5.*sim->U[i][j+2] + 15.*sim->U[i][j+1] - 16 * 0.);

        j = sim->ny+1;
        sim->U[i][j] = -0.2 * (sim->U[i][j-3] - 5.*sim->U[i][j-2] + 15.*sim->U[i][j-1] - 16 * 0.);
    }
#   endif

#   if TEMPERATURE
    double T_wall;
    for (i = 1; i < sim->nx + 1; i++) {
        
        j = 0;
        T_wall =  0.5 * (1. - tanh( ((i - 0.5) * sim->h - L_/2.) / 1.)) * TMAX;
        sim->T[i][j] = -0.2 * (sim->T[i][j+3] - 5.*sim->T[i][j+2] + 15.*sim->T[i][j+1] - 16.*T_wall);

        j = sim->ny+1;
        T_wall = 0.5 * (1. - tanh( ((i - 0.5) * sim->h - L_/2.) / 1.)) * TMIN;
        sim->T[i][j] = -0.2 * (sim->T[i][j-3] - 5.*sim->T[i][j-2] + 15.*sim->T[i][j-1] - 16.*T_wall);
    }
#   endif
}


void save_fields(Sim_data *sim, int t) {
    FILE *ptr, *ptr_u, *ptr_v, *ptr_p, *ptr_T;
    int i;

    if (t == 0) {
        ptr = fopen(filename_params, "w");
        fprintf(ptr, "%d %d %d %d %d\n", sim->nt / SAVE_MODULO, sim->nx, sim->ny, sim->n, SAVE_MODULO);
        fprintf(ptr, "%lf %lf %lf %lf %d %d %d %d %d\n", RE, sim->tsim, sim->dt, sim->h, L_, H_, LBOX, D_IN, D_BOT);
        fprintf(ptr, "%lf %lf %lf %lf %lf %lf %d\n", ALPHA, STROUHAL, SIWNG_START, KAPPA_Y, STROUHAL_Y, PERT_START, N_CYCLES);
        fprintf(ptr, "%d %lf %lf %lf %lf\n", TEMPERATURE, PR, GR, TMAX, TMIN);
        fclose(ptr);
    }
    
    if (t == 0) {
        ptr_u = fopen(filename_u, "w");
        ptr_v = fopen(filename_v, "w");
        ptr_p = fopen(filename_p, "w");
        ptr_T = fopen(filename_T, "w");
    } else {
        ptr_u = fopen(filename_u, "a");
        ptr_v = fopen(filename_v, "a");
        ptr_p = fopen(filename_p, "a");
        ptr_T = fopen(filename_T, "a");
    }

    for (i = 0; i < sim->size_u; i++) fprintf(ptr_u, FMT, sim->u_data[i]);  // velocity x
    for (i = 0; i < sim->size_v; i++) fprintf(ptr_v, FMT, sim->v_data[i]);  // velocity y
    for (i = 0; i < sim->size_p; i++) fprintf(ptr_p, FMT, sim->p_data[i]);  // pressure

#   if TEMPERATURE
    int j;
    for (i = 1; i < sim->nx + 1; i++) {
        for (j = 1; j < sim->ny + 1; j++) {
            fprintf(ptr_T, FMT, sim->T[i][j]);  // temperature
        }
    }
#   endif

    fclose(ptr_u);
    fclose(ptr_v);
    fclose(ptr_p);
    fclose(ptr_T);
 }


void set_bd_conditions(Sim_data *sim, double **U, double **V) {
    int i, j;
    double coef;
    
    /*// Inflow boundary condition for u  // never changes
    for (j = 1; j <= sim->ny; j++) {   // Inflow uniform profile u
        U[0][j] = 1.;   // could do a special profile
    }*/

    // Outflow boundary condition for u
    i = sim->nx;
    coef = sim->dt / sim->h * (1. - sim->uMesh);
    for (j = 1; j <= sim->ny; j++) {  // u is advected at velocity (1 - uMesh)
        U[i][j] -= coef * (sim->U[i][j] - sim->U[i-1][j]);
    }

#   if NO_SLIP
    // Channel walls: condition for v
    for (i = 0; i < sim->nx + 2; i++) { // no-through flow
        V[i][0] = sim->vMesh;             // below
        V[i][sim->ny] = sim->vMesh;       // above
    }
#   endif


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
    }
}

void set_ghost_points(Sim_data *sim) {
    int i, j;
    double coef, w_last, w_left;

    // External walls
#   if NO_SLIP  // (v = 0 and u = 0)
    for (i = 0; i < sim->nx + 1; i++) {
        j = 0;
        sim->U[i][j] = -0.2 * (sim->U[i][j+3] - 5.*sim->U[i][j+2] + 15.*sim->U[i][j+1] - 16.*sim->uMesh);
        j = sim->ny+1;
        sim->U[i][j] = -0.2 * (sim->U[i][j-3] - 5.*sim->U[i][j-2] + 15.*sim->U[i][j-1] - 16.*sim->uMesh);
    }
#   else  // (v = 0 and w = 0 --> du/dy = 0)
    for (i = 0; i < sim->nx + 1; i++) {                // zero voriticty
        sim->U[i][0] = sim->U[i][1];                   // below
        sim->U[i][sim->ny + 1] = sim->U[i][sim->ny];   // above
    }
#   endif

    // Inflow boundary condition for v
    for (j = 1; j < sim->ny; j++) {  // trick with ghost for zero v at inflow
        sim->V[0][j] = -0.2 * (sim->V[3][j] - 5. * sim->V[2][j] + 15. * sim->V[1][j] - 16. * 0.);  // 0 or vMesh ?
    }
    
    // Outflow boundary condition for v
    i = sim->nx;
    coef = sim->dt / sim->h * (1. - sim->uMesh);
    for (j = 1; j < sim->ny; j++) {
        w_last = (sim->V[i+1][j] - sim->V[i  ][j]) - (sim->U[i  ][j+1] - sim->U[i  ][j]);
        w_left = (sim->V[i  ][j] - sim->V[i-1][j]) - (sim->U[i-1][j+1] - sim->U[i-1][j]);
        sim->V[i+1][j] = (sim->V[i][j] + sim->U[i][j+1] - sim->U[i][j]) + w_last - coef * (w_last - w_left);
        // [ v_(i+0.5,j) - v_(i-0.5,j) - u_(i,j+0.5) + u(i,j-0.5) ] / h = w_ij + dt uc / h * [ w_ij - w_(i-1, j) ]
    }


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
        sim->U[i][j] = -0.2 * (sim->U[i][j-3] - 5.*sim->U[i][j-2] + 15.*sim->U[i][j-1] - 16.*sim->uMesh);
        j = (D_BOT + 1) * sim->n;  // above
        sim->U[i][j] = -0.2 * (sim->U[i][j+3] - 5.*sim->U[i][j+2] + 15.*sim->U[i][j+1] - 16.*sim->uMesh);
    }
}

void set_boundary_temperature(Sim_data *sim) {
    int i, j;
    double coef, T_wall;
    coef = sim->dt / sim->h * (1. - sim->uMesh);

    for (i = 1; i < sim->nx + 1; i++) {  // Above and below (DIRICHLET)        
        j = 0;
        T_wall = 0.5 * (1. - tanh( ((i - 0.5) * sim->h - L_/2.) / 1.)) * TMAX;
        sim->T[i][j] = -0.2 * (sim->T[i][j+3] - 5.*sim->T[i][j+2] + 15.*sim->T[i][j+1] - 16. * T_wall);

        j = sim->ny+1;
        T_wall = 0.5 * (1. - tanh( ((i - 0.5) * sim->h - L_/2.) / 1.)) * TMIN;
        sim->T[i][j] = -0.2 * (sim->T[i][j-3] - 5.*sim->T[i][j-2] + 15.*sim->T[i][j-1] - 16. * T_wall);
    }

    for (j = 1; j < sim->ny + 1; j++) sim->T[0        ][j] = sim->T[1][j];                                         // inflow
    for (j = 1; j < sim->ny + 1; j++) sim->T[sim->nx+1][j] -= coef * (sim->T[sim->nx+1][j] - sim->T[sim->nx][j]);  // outflow (advected)

    i = sim->n * D_IN + 1;  // left wall of rectanle
    for (j = sim->n * D_BOT + 1; j < sim->n * (D_BOT + 1) + 1; j++)
        sim->T[i][j] = sim->T[i-1][j];
    
    i = sim->n * (D_IN + LBOX);  // right wall of rectanle
    for (j = sim->n * D_BOT + 1; j < sim->n * (D_BOT + 1) + 1; j++)
        sim->T[i][j] = sim->T[i+1][j];

    // ! Ghost points in the corners are overwritten
    for (i = sim->n * D_IN + 1; i < sim->n * (D_IN + LBOX) + 1; i++) {
        j = sim->n * D_BOT + 1;  // lower wall of rectangle
        sim->T[i][j] = sim->T[i][j-1];
        j = sim->n * (D_BOT + 1);  // upper wall of rectangle
        sim->T[i][j] = sim->T[i][j+1];
    }
}


void compute_convection(Sim_data *sim) {
    int i, j;

    double **U = sim->U;
    double **V = sim->V;
    double uMesh = sim->uMesh;
    double vMesh = sim->vMesh;

    double factor = 1. / (4. * sim->h);

#   if (CONVECTION_MODE == 2)
    factor *= 0.5;
#   endif

    double **H = sim->HX;

    int i_s, i_f, j_s, j_f;

    for (int block = 0; block < 4; block++) {
        i_s = sim->i_start[block];
        i_f = sim->i_final[block];
        j_s = sim->j_start[block];
        j_f = sim->j_final[block];

        for (i = i_s;  i < i_f; i++) {
            for (j = j_s; j < j_f; j++) {

                H[i][j] = factor * (   // equivalent to divide by 4*h
#               if (CONVECTION_MODE == 0) || (CONVECTION_MODE == 2)  // Advective form for u field
                        + (U[i+1][j  ] + U[i  ][j  ] - 2.*uMesh) * (U[i+1][j  ] - U[i  ][j  ])
                        + (U[i  ][j  ] + U[i-1][j  ] - 2.*uMesh) * (U[i  ][j  ] - U[i-1][j  ])
                        + (V[i+1][j  ] + V[i  ][j  ] - 2.*vMesh) * (U[i  ][j+1] - U[i  ][j  ])
                        + (V[i+1][j-1] + V[i  ][j-1] - 2.*vMesh) * (U[i  ][j  ] - U[i  ][j-1])
#               endif
#               if  (CONVECTION_MODE == 1) || (CONVECTION_MODE == 2)  // Divergence form for u field
                        + (U[i+1][j  ] + U[i  ][j  ] - 2.*uMesh) * (U[i+1][j  ] + U[i  ][j  ])
                        - (U[i  ][j  ] + U[i-1][j  ] - 2.*uMesh) * (U[i  ][j  ] + U[i-1][j  ])
                        + (V[i+1][j  ] + V[i  ][j  ] - 2.*vMesh) * (U[i  ][j+1] + U[i  ][j  ])
                        - (V[i+1][j-1] + V[i  ][j-1] - 2.*vMesh) * (U[i  ][j  ] + U[i  ][j-1])
#               endif
                );
            }
        }
    }

    H = sim->HY;
    for (int block = 4; block < 8; block++) {
        i_s = sim->i_start[block];
        i_f = sim->i_final[block];
        j_s = sim->j_start[block];
        j_f = sim->j_final[block];

        for (i = i_s;  i < i_f; i++) {
            for (j = j_s; j < j_f; j++) {
                
                H[i][j] = factor * (
#               if (CONVECTION_MODE == 0) || (CONVECTION_MODE == 2)  // Advective form for v field
                        + (U[i  ][j+1] + U[i  ][j  ] - 2.*uMesh) * (V[i+1][j  ] - V[i  ][j  ])
                        + (U[i-1][j+1] + U[i-1][j  ] - 2.*uMesh) * (V[i  ][j  ] - V[i-1][j  ])
                        + (V[i  ][j+1] + V[i  ][j  ] - 2.*vMesh) * (V[i  ][j+1] - V[i  ][j  ])
                        + (V[i  ][j  ] + V[i  ][j-1] - 2.*vMesh) * (V[i  ][j  ] - V[i  ][j-1])
#               endif
#               if  (CONVECTION_MODE == 1) || (CONVECTION_MODE == 2)  // Divergence form for v field
                        + (U[i  ][j+1] + U[i  ][j  ] - 2.*uMesh) * (V[i+1][j  ] + V[i  ][j  ])
                        - (U[i-1][j+1] + U[i-1][j  ] - 2.*uMesh) * (V[i  ][j  ] + V[i-1][j  ])
                        + (V[i  ][j+1] + V[i  ][j  ] - 2.*vMesh) * (V[i  ][j+1] + V[i  ][j  ])
                        - (V[i  ][j  ] + V[i  ][j-1] - 2.*vMesh) * (V[i  ][j  ] + V[i  ][j-1])
#               endif
                );

            }
        }
    }
}

void compute_convection_temperature(Sim_data *sim) {
    int i, j;

    double **T = sim->T;
    double **HT = sim->HT;
    double **U = sim->U;
    double **V = sim->V;

    double uMesh = sim->uMesh;
    double vMesh = sim->vMesh;

    double factor = 1. / (2. * sim->h);

#   if (CONVECTION_MODE == 2)
    factor *= 0.5;
#   endif

    int i_s, i_f, j_s, j_f;

    for (int block = 8; block < 12; block++) {
        i_s = sim->i_start[block];
        i_f = sim->i_final[block];
        j_s = sim->j_start[block];
        j_f = sim->j_final[block];

        for (i = i_s;  i < i_f; i++) {
            for (j = j_s; j < j_f; j++) {

                HT[i][j] = factor * (
#               if (CONVECTION_MODE == 0) || (CONVECTION_MODE == 2)  // Advective form for v field
                        + (U[i  ][j  ] - uMesh) * (T[i+1][j  ] - T[i  ][j  ])
                        + (U[i-1][j  ] - uMesh) * (T[i  ][j  ] - T[i-1][j  ])
                        + (V[i  ][j  ] - vMesh) * (T[i  ][j+1] - T[i  ][j  ])
                        + (V[i  ][j-1] - vMesh) * (T[i  ][j  ] - T[i  ][j-1])
#               endif
#               if  (CONVECTION_MODE == 1) || (CONVECTION_MODE == 2)  // Divergence form for v field
                        + (U[i  ][j  ] - uMesh) * (T[i+1][j  ] + T[i  ][j  ])
                        - (U[i-1][j  ] - uMesh) * (T[i  ][j  ] + T[i-1][j  ])
                        + (V[i  ][j  ] - vMesh) * (T[i  ][j+1] + T[i  ][j  ])
                        - (V[i  ][j-1] - vMesh) * (T[i  ][j  ] + T[i  ][j-1])
#               endif
                );

            }
        }
    }
}


void predictor_step(Sim_data *sim) {
    int i, j;
    // double conv, pres, diff;

    double **U = sim->U;
    double **V = sim->V;

    double coef_1 = (sim->t == 0) ? -1. : -1.5;
    double coef_2 = (sim->t == 0) ? +0. : +0.5;
    double alpha = 1. / (RE * sim->h * sim->h);
#   if TEMPERATURE
    double beta = 0.5 * GR / (RE * RE);
#   endif

    int i_s, i_f, j_s, j_f;

    // update u field
    for (int block = 0; block < 4; block++) {
        i_s = sim->i_start[block];
        i_f = sim->i_final[block];
        j_s = sim->j_start[block];
        j_f = sim->j_final[block];

        for (i = i_s;  i < i_f; i++) {
            for (j = j_s; j < j_f; j++) {
                            
                // conv = coef_1 * sim->HX[i][j] + coef_2 * sim->HX_[i][j];  
                // pres = -(sim->P[i][j-1] - sim->P[i-1][j-1]) * sim->n; 
                // diff = alpha * (U[i+1][j] + U[i-1][j] - 4*U[i][j] + U[i][j+1] + U[i][j-1]);  

                sim->US[i][j] = U[i][j] + sim->dt * (
                    + coef_1 * sim->HX[i][j] + coef_2 * sim->HX_[i][j]                     // (u du/dx + v du/dy) at t_n and t_(n-1)
                    - (sim->P[i][j-1] - sim->P[i-1][j-1]) * sim->n                         // dp/dx
                    + alpha * (U[i+1][j] + U[i-1][j] - 4*U[i][j] + U[i][j+1] + U[i][j-1])  // d2u/dx2 + d2u/dy2
                );
            }
        }
    }

    // update v field
    // double count1 = 0.; double count2 = 0.; double count3 = 0.; double count4 = 0.; double c5 = 0.;
    for (int block = 4; block < 8; block++) {
        i_s = sim->i_start[block];
        i_f = sim->i_final[block];
        j_s = sim->j_start[block];
        j_f = sim->j_final[block];

        for (i = i_s;  i < i_f; i++) {
            for (j = j_s; j < j_f; j++) {
                
                // conv = coef_1 * sim->HY[i][j] + coef_2 * sim->HY_[i][j];  
                // pres = -(sim->P[i-1][j] - sim->P[i-1][j-1]) * sim->n;  
                // diff = alpha * (V[i+1][j] + V[i-1][j] - 4*V[i][j] + V[i][j+1] + V[i][j-1]);  

                sim->VS[i][j] = V[i][j] + sim->dt * (
                    + coef_1 * sim->HY[i][j] + coef_2 * sim->HY_[i][j]                    // (u dv/dx + v dv/dy) at t_n and t_(n-1)
                    - (sim->P[i-1][j] - sim->P[i-1][j-1]) * sim->n                        // dp/dy
                    + alpha * (V[i+1][j] + V[i-1][j] - 4*V[i][j] + V[i][j+1] + V[i][j-1]) // d2v/dx2 + d2v/dy2
#                   if TEMPERATURE
                    + beta * (sim->T[i][j+1] + sim->T[i][j])  // AVG(T above, T below)
#                   endif
                );

                // count1 += coef_1 * sim->HY[i][j] + coef_2 * sim->HY_[i][j];
                // count2 += (sim->P[i-1][j] - sim->P[i-1][j-1]) * sim->n;
                // count3 += alpha * (V[i+1][j] + V[i-1][j] - 4*V[i][j] + V[i][j+1] + V[i][j-1]);
                // count4 += beta * (sim->T[i][j+1] + sim->T[i][j]);
                // c5 += 1.;
            }
        }
    }
    // printf("Count 1 = %7.3le  Count 2 = %7.3le  Count 3 = %7.3le  Count 4 = %7.3le\n", count1/c5, count2/c5, count3/c5, count4/c5);
}


void corrector_step(Sim_data *sim) {
    int i, j;
    double dphi;

    int i_s, i_f, j_s, j_f;

    // update u field
    for (int block = 0; block < 4; block++) {
        i_s = sim->i_start[block];
        i_f = sim->i_final[block];
        j_s = sim->j_start[block];
        j_f = sim->j_final[block];
        
        for (i = i_s;  i < i_f; i++) {
            for (j = j_s; j < j_f; j++) {
                dphi = (sim->PHI[i][j-1] - sim->PHI[i-1][j-1]) * sim->n;  // equivalent to divide by h
                sim->U[i][j] = sim->US[i][j] - sim->dt * dphi;  // - d phi / dx
            }
        }
    }

    // update v field
    for (int block = 4; block < 8; block++) {
        i_s = sim->i_start[block];
        i_f = sim->i_final[block];
        j_s = sim->j_start[block];
        j_f = sim->j_final[block];
        
        for (i = i_s;  i < i_f; i++) {
            for (j = j_s; j < j_f; j++) {
                dphi = (sim->PHI[i-1][j] - sim->PHI[i-1][j-1]) * sim->n;
                sim->V[i][j] = sim->VS[i][j] - sim->dt * dphi; // - d phi / dy
            }
        }
    }

    // update p field
    for (i = 0, j = sim->size_p; i < sim->size_p; i++, j++) {
        sim->p_data[i] += sim->p_data[j];  // p = p + phi
    }

    // update T field
#   if TEMPERATURE
    double coef_1 = (sim->t == 0) ? -1. : -1.5;
    double coef_2 = (sim->t == 0) ? +0. : +0.5;
    double diff = 1./(PR * RE) * 1. / (sim->h * sim->h);  // 1/(PR RE) * 1/h^2

    // compute laplacian of temperature
    // use existing array PHI for this purpose -> index -1 since size of PHI smaller
    for (i = 1; i < sim->nx + 1; i++) {
        for (j = 1; j < sim->ny + 1; j++) {
            sim->PHI[i-1][j-1] = sim->T[i+1][j] + sim->T[i-1][j] + sim->T[i][j+1] + sim->T[i][j-1] - 4. * sim->T[i][j];
        }
    }
    // ! Handle 4 points where ghost points were overwritten (corners of the rectangle)
    i = (D_IN) * sim->n;
    j = (D_BOT) * sim->n + 1;
    sim->PHI[i-1][j-1] = sim->T[i-1][j] + sim->T[i][j+1] + sim->T[i][j-1] - 3. * sim->T[i][j]; // T[i+1] = T[i]
    j = (D_BOT + 1) * sim->n;
    sim->PHI[i-1][j-1] = sim->T[i-1][j] + sim->T[i][j+1] + sim->T[i][j-1] - 3. * sim->T[i][j]; // T[i+1] = T[i]
    
    i = (D_IN + LBOX) * sim->n + 1;
    j = (D_BOT) * sim->n + 1;
    sim->PHI[i-1][j-1] = sim->T[i+1][j] + sim->T[i][j+1] + sim->T[i][j-1] - 3. * sim->T[i][j]; // T[i-1] = T[i]
    j = (D_BOT + 1) * sim->n;
    sim->PHI[i-1][j-1] = sim->T[i+1][j] + sim->T[i][j+1] + sim->T[i][j-1] - 3. * sim->T[i][j]; // T[i-1] = T[i]
    
    for (i = 1; i < sim->nx + 1; i++) {
        for (j = 1; j < sim->ny + 1; j++) {
            sim->T[i][j] += sim->dt * (coef_1 * sim->HT[i][j] + coef_2 * sim->HT_[i][j] + diff * sim->PHI[i-1][j-1]);
        }
    }
#   endif
}


void swap_next_previous(Sim_data *sim) {
    double **tmp;
    
    tmp = sim->HX;
    sim->HX = sim->HX_;
    sim->HX_ = tmp;

    tmp = sim->HY;
    sim->HY = sim->HY_;
    sim->HY_ = tmp;

#   if TEMPERATURE
    tmp = sim->HT;
    sim->HT = sim->HT_;
    sim->HT_ = tmp;
#   endif
}


void set_mesh_velocity(Sim_data *sim) {
    double t_now = sim->dt * sim->t;
    
    if (SIWNG_START < t_now) {
        sim->uMesh = ALPHA * sin(2. * M_PI * STROUHAL * (t_now - SIWNG_START));
    } else {
        sim->uMesh = 0.;
    }

    if ((PERT_START <= t_now) && (t_now <= PERT_START + ((double) N_CYCLES / STROUHAL_Y))) {
        sim->vMesh = KAPPA_Y * 2. * M_PI * STROUHAL_Y * sin(2. * M_PI * STROUHAL_Y * (t_now - PERT_START));
    } else {
        sim->vMesh = 0.;
    }
}


void check_boundary(Sim_data *sim) {
    int i, j;

    // set identity for upper and lower walls
    for (i = 0; i < sim->nx + 1; i++) {
        for (j = 0; j < sim->ny + 2; j += sim->ny+1) {
            printf("u[%d][%d] = %lf \t us = %lf\n", i, j, sim->U[i][j], sim->US[i][j]);
        }
    }

    // set identity for left and right walls
    for (j = 1; j < sim->ny + 1; j++) {
        for (i = 0; i < sim->nx+1; i += sim->nx) {
            printf("u[%d][%d] = %lf \t us = %lf\n", i, j, sim->U[i][j], sim->US[i][j]);
        }
    }

    // set identity for rectangle
    for (i = sim->n * D_IN; i < (D_IN + LBOX) * sim->n + 1; i++) {
        for (j = D_BOT * sim->n + 1; j < (D_BOT + 1) * sim->n + 1; j++) {
            printf("u[%d][%d] = %lf \t us = %lf\n", i, j, sim->U[i][j], sim->US[i][j]);
        }
    }

    // set identity for upper and lower walls
    for (i = 0; i < sim->nx + 2; i++) {
        for (j = 0; j < sim->ny + 1; j += sim->ny) {
            printf("v[%d][%d] = %lf \t vs = %lf\n", i, j, sim->V[i][j], sim->VS[i][j]);
        }
    }

    // set identity for left and right walls
    for (j = 1; j < sim->ny; j++) {
        for (i = 0; i < sim->nx+2; i += sim->nx + 1) {
            printf("v[%d][%d] = %lf \t vs = %lf\n", i, j, sim->V[i][j], sim->VS[i][j]);
        }
    }

    // set identity for rectangle
    for (i = sim->n * D_IN + 1; i < (D_IN + LBOX) * sim->n + 1; i++) {
        for (j = D_BOT * sim->n; j < (D_BOT + 1) * sim->n + 1; j++) {
            printf("v[%d][%d] = %lf \t vs = %lf\n", i, j, sim->V[i][j], sim->VS[i][j]);
        }
    }
}


void free_Sim_data(Sim_data *sim) {
    free(sim->u_data); free(sim->U);
    free(sim->v_data); free(sim->V);
    free(sim->p_data); free(sim->P);
#   if TEMPERATURE
    free(sim->T_data); free(sim->T);
#   endif
    free(sim);
}
