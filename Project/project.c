#include "poisson.h"
#include "project.h"

char *myPath = "./data/";
char filename_params[50];
char filename_u[50];
char filename_v[50];
char filename_p[50];
char filename_w[50];


void init_data_sim(data_Sim *sim) {
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

#if DEBUG
    sim->n = 4;
#else
    sim->n = 10;  // 1. / h
#endif

    sim->nt = 100;
    sim->nx = L_ * sim->n;
    sim->ny = H_ * sim->n;

    // sim->dt = TEND / (double) sim->nt;
    sim->dt = 0.01;
    sim->h = H_ / (double) sim->ny;
    
    sim->uMesh = 0.;
    sim->vMesh = 0.;

    int size_u = (sim->nx + 1) * (sim->ny + 2);
    int size_v = (sim->nx + 2) * (sim->ny + 1);
    int size_p = (sim->nx) * (sim->ny);
    sim->size_u = size_u;
    sim->size_v = size_v;
    sim->size_p = size_p;
    
    sim->u = (double *)calloc(5*size_u, sizeof(double));
    sim->u_star = sim->u + size_u;
    //sim->u_prev = sim->u + 2*size_u;
    sim->Hx = sim->u + 3*size_u;
    sim->Hx_prev = sim->u + 4*size_u;

    sim->v = (double *)calloc(5*size_v, sizeof(double));
    sim->v_star = sim->v + size_v;
    //sim->v_prev = sim->v + 2*size_v;
    sim->Hy = sim->v + 3*size_v;
    sim->Hy_prev = sim->v + 4*size_v;

    sim->p = (double *)calloc(2*size_p, sizeof(double));
    sim->phi = sim->p + size_p;
}



void init_fields(data_Sim *sim) {
    int i, j;
    int k = sim->ny + 2;

    int i_w_left =  D_IN * sim->n;           // on boundary
    int i_w_right = (D_IN + LBOX) * sim->n;  // on boundary
    int j_w_below = D_BOT * sim->n + 1;      // in rectangle
    int j_w_above = (D_BOT + 1) * sim->n;    // in rectangle


    // set u to Uinf at beginning (0 is good for v)
    for (i = 0; i < sim->size_u; i++) {
        sim->u[i] = 1.;
        sim->u_star[i] = 1.;
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
        sim->u[i*k + j] = -0.2 * 11.; // since u(t=0) is 1.

        j = j_w_above;
        sim->u[i*k + j] = -0.2 * 11.; // since u(t=0) is 1.
    }
}


void save_fields(data_Sim *sim, int t) {
    FILE *ptr, *ptr_u, *ptr_v, *ptr_p, *ptr_w;

    if (t == 0) {
        ptr = fopen(filename_params, "w");
        fprintf(ptr, "%d %d %d\n", sim->nt, sim->nx, sim->ny);
        fprintf(ptr, "%lf %d %d %d %d %d %lf %lf %lf %lf\n", TEND, L_, H_, LBOX, D_IN, D_BOT, sim->dt, sim->h, ALPHA, STROUHAL);
        fclose(ptr);
    }

    int i, j;
    int k = sim->ny+2;
    int l = sim->ny+1;
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
            
            u_up = sim->u[i*k + j + 1]; // u above
            u_dw = sim->u[i*k + j];     // u below
            v_lf = sim->v[i*l + j];     // v left
            v_rg = sim->v[(i+1)*l + j]; // v right
            
            fprintf(ptr_u, FMT, 0.5 * (u_up + u_dw));                       // u average 
            fprintf(ptr_v, FMT, 0.5 * (v_lf + v_rg));                       // v average
            fprintf(ptr_w, FMT, ((v_rg - v_lf) - (u_up - u_dw)) / sim->h);  // vorticity
        }
    }

    // loop over centers
    for (i = 0; i < sim->size_p; i++) {
        fprintf(ptr_p, FMT, sim->p[i]);  // pressure
    }
    
    fprintf(ptr_u, "\n");    fclose(ptr_u);
    fprintf(ptr_v, "\n");    fclose(ptr_v);
    fprintf(ptr_w, "\n");    fclose(ptr_w);    
    fprintf(ptr_p, "\n");    fclose(ptr_p);

}


void set_bd_conditions(data_Sim *sim, double *u, double *v) {
    int i, j;
    int idx_u;
    
    int k = sim->ny + 2;
    int l = sim->ny + 1;

    /*// Inflow boundary condition for u
    for (j = 1; j <= sim->ny; j++) {   // Inflow uniform profile u
        u[j] = 1.;   // could do a special profile
    }*/

    // Outflow boundary condition for u
    for (j = 1; j <= sim->ny; j++) {  // u is advected at velocity (1 - uMesh)
        idx_u = sim->nx * k + j;
        u[idx_u] += sim->dt * (1. - sim->uMesh) * (sim->u[idx_u] - LL(sim->u, idx_u, k)) / sim->h;
    }

    /*// Lateral walls condition for v
    for (i = 0; i < sim->nx + 2; i++) { // no-through flow
        v[i*l] = 0.;               // below
        v[i*l + (l-1)] = 0.;       // above
    }*/

    // Left wall of rectangle (u)
    i = D_IN * sim->n;
    for (j = D_BOT * sim->n + 1; j <= (D_BOT + 1) * sim->n; j++)
        u[i*k + j] = sim->uMesh;
    
    // Right wall of rectangle (u)
    i = (D_IN + LBOX) * sim->n;
    for (j = D_BOT * sim->n + 1; j <= (D_BOT + 1) * sim->n; j++)
        u[i*k + j] = sim->uMesh;

    // Upper and lower walls of rectangle (v)
    for (i = D_IN * sim->n + 1; i < (D_IN + LBOX) * sim->n; i++) {
        
        // lower wall of the rectangle
        j = D_BOT * sim->n;
        v[i*l + j] = sim->vMesh;

        // upper wall of the rectangle
        j = (D_BOT + 1) * sim->n;
        v[i*l + j] = sim->vMesh;
    }
}

void set_ghost_points(data_Sim *sim) {
    int i, j;
    int idx_u, idx_v;
    
    int k = sim->ny + 2;
    int l = sim->ny + 1;
    
    // External walls (v = 0 and w = 0 --> du/dy = 0)
    for (i = 0; i < sim->nx + 1; i++) {                // zero voriticty
        sim->u[i*k] = sim->u[i*k + 1];                 // below
        sim->u[i*k + (k-1)] = sim->u[i*k + (k-2)];     // above
    }

    // Inflow boundary condition for v
    for (j = 1; j < sim->ny; j++) {  // trick with ghost for zero v at inflow
        sim->v[j] = -0.2 * (sim->v[3*l + j] - 5. * sim->v[2*l + j] + 15. * sim->v[l + j] - 16. * 0.);  // 0 or vMesh ?
    }
    
    // Outflow boundary condition for v
    double w_last, w_left;
    for (j = 1; j < sim->ny; j++) {
        idx_v = sim->nx * l + j;
        idx_u = sim->nx * k + j;
        w_last = sim->v[idx_v + l] - sim->v[idx_v] - sim->u[idx_u+1] + sim->u[idx_u];
        w_left = sim->v[idx_v] - sim->v[idx_v-l] - sim->u[idx_u-k+1] + sim->u[idx_u-k];
        
        sim->v[idx_v + l] = (sim->v[idx_v] + sim->u[idx_u+1] - sim->u[idx_u]) \
                           + w_last + sim->dt / sim->h * (1 - sim->uMesh) * (w_last - w_left);

        // [ v_(i+0.5,j) - v_(i-0.5,j) - u_(i,j+0.5) + u(i,j-0.5) ] / h = w_ij + dt uc / h * [ w_ij - w_(i-1, j) ]
    }

    // Left wall of rectangle
    i = D_IN * sim->n + 1;
    for (j = D_BOT * sim->n + 1; j < (D_BOT + 1) * sim->n; j++)
        sim->v[i*l + j] = -0.2 * (sim->v[(i-3)*l + j] - 5.*sim->v[(i-2)*l + j] + 15.*sim->v[(i-1)*l + j] - 16.*sim->vMesh);

    // Right wall of rectangle
    i = (D_IN + LBOX) * sim->n;
    for (j = D_BOT * sim->n + 1; j < (D_BOT + 1) * sim->n; j++)
        sim->v[i*l + j] = -0.2 * (sim->v[(i+3)*l + j] - 5.*sim->v[(i+2)*l + j] + 15.*sim->v[(i+1)*l + j] - 16.*sim->vMesh);

    // Upper and lower walls of rectangle
    for (i = D_IN * sim->n + 1; i < (D_IN + LBOX) * sim->n; i++) {
    
        j = D_BOT * sim->n + 1;  // below
        sim->u[i*k + j] = -0.2 * (sim->u[i*k + j - 3] - 5.*sim->u[i*k + j - 2] + 15 * sim->u[i*k + j - 1] - 16*sim->uMesh);

        j = (D_BOT + 1) * sim->n;  // above
        sim->u[i*k + j] = -0.2 * (sim->u[i*k + j - 3] - 5.*sim->u[i*k + j - 2] + 15 * sim->u[i*k + j - 1] - 16*sim->uMesh);
    }
}


void compute_convection(data_Sim *sim) {
    int i, j, k, l, idx_u, idx_v;
    k = sim->ny + 2;
    l = sim->ny + 1;

    double *u = sim->u;
    double *v = sim->v;
    double uMesh = sim->uMesh;
    double vMesh = sim->vMesh;

    double *H = sim->Hx;
    for (i = 1; i < sim->nx; i++) {
        for (j = 1; j <= sim->ny; j++) {
            idx_u = i * k + j;
            idx_v = i * l + j; // up left wrt u

            // Advective form for u field
            H[idx_u] = (RR(u, idx_u, k) + u[idx_u] - 2.*uMesh) * (RR(u, idx_u, k) - u[idx_u]) \
                     + (u[idx_u] + LL(u, idx_u, k) - 2.*uMesh) * (u[idx_u] - LL(u, idx_u, k)) \
                     + (AL(v, idx_v, l) + AR(v, idx_v, l) - 2.*vMesh) * (AA(u, idx_u, k) - u[idx_u]) \
                     + (BL(v, idx_v, l) + BR(v, idx_v, l) - 2.*vMesh) * (u[idx_u] - BB(u, idx_u, k));
            
            // Divergence form for u field
            H[idx_u] += (RR(u, idx_u, k) + u[idx_u]) * (RR(u, idx_u, k) + u[idx_u] - 2.*uMesh) \
                      - (LL(u, idx_u, k) + u[idx_u]) * (LL(u, idx_u, k) + u[idx_u] - 2.*uMesh) \
                      + (AA(u, idx_u, k) + u[idx_u]) * (AL(v, idx_v, l) + AR(v, idx_v, l) - 2.*vMesh) \
                      - (BB(u, idx_u, k) + u[idx_u]) * (BL(v, idx_v, l) + BR(v, idx_v, l) - 2.*vMesh);

            H[idx_u] /= 2. * (4. * sim->h);  // possible since dx = dy = h; factor "2" since avg of Adv / Div
        }
    }

    H = sim->Hy;
    for (i = 1; i <= sim->nx; i++) {
        for (j = 1; j < sim->ny; j++) {
            idx_u = (i-1) * k + (j + 1);  // up left wrt v
            idx_v = i * l + j; 
            
            // Advective form for v field
            H[idx_v] = (AR(u, idx_u, k) + BR(u, idx_u, k) - 2.*uMesh) * (RR(v, idx_v, l) - v[idx_v]) \
                     + (AL(u, idx_u, k) + BL(u, idx_u, k) - 2.*uMesh) * (v[idx_v] - LL(v, idx_v, l)) \
                     + (AA(v, idx_v, l) + v[idx_v] - 2.*vMesh) * (AA(v, idx_v, l) - v[idx_v]) \
                     + (BB(v, idx_v, l) + v[idx_v] - 2.*vMesh) * (v[idx_v] - BB(v, idx_v, l));
            
            // Divergence form for v field
            H[idx_v] += (RR(v, idx_v, l) + v[idx_v]) * (AR(u, idx_u, k) + BR(u, idx_u, k) - 2.*uMesh) \
                      - (LL(v, idx_v, l) + v[idx_v]) * (AL(u, idx_u, k) + BL(u, idx_u, k) - 2.*uMesh) \
                      + (AA(v, idx_v, l) + v[idx_v]) * (AA(v, idx_v, l) + v[idx_v] - 2.*vMesh) \
                      - (BB(v, idx_v, l) + v[idx_v]) * (BB(v, idx_v, l) + v[idx_v] - 2.*vMesh);

            H[idx_v] /= 2. * (4. * sim->h);
        }
    }
}


void save_debug(data_Sim *sim) {
    FILE *ptr_u = fopen("fields_u.txt", "w");
    FILE *ptr_v = fopen("fields_v.txt", "w");
    
    int i, j;
    int k = sim->ny + 2;
    for (j = sim->ny+1; j >= 0; j--) {
        for (i = 0 ; i < sim->nx+1; i++) {
            fprintf(ptr_u, "%3.0lf/%7.2lf     ", sim->u[i*k+j], sim->Hx[i*k+j]);
        }
        fprintf(ptr_u, "\n");
    }

    k = sim->ny + 1;
    for (j = sim->ny; j >= 0; j--) {
        for (i = 0 ; i <= sim->nx+1; i++) {
            fprintf(ptr_v, "%3.0lf/%7.2lf     ", sim->v[i*k+j], sim->Hy[i*k+j]);
        }
        fprintf(ptr_v, "\n");
    }

    fclose(ptr_u);
    fclose(ptr_v);
    printf("done\n");
}


void update_1(data_Sim *sim, int t) {
    int i, j, idx_u, idx_v, idx_p;
    double conv, pres, diff;
    
    int k = sim->ny + 2;
    int l = sim->ny + 1;

    double *u = sim->u;
    double *v = sim->v;

    double coef_1 = (t = 1) ? -1. : -1.5;
    double coef_2 = (t = 1) ? +0. : +0.5;
    double alpha = 1. / (RE * sim->h * sim->h);

    // update u field
    for (i = 1; i < sim->nx; i++) {
        for (j = 1; j <= sim->ny; j++) {
            
            idx_u = i * k + j;
            idx_p = i * sim->ny + j-1;
            
            conv = coef_1 * sim->Hx[idx_u] + coef_2 * sim->Hx_prev[idx_u];  // (u du/dx + v du/dy) at t_n and t_(n-1)
            pres = -(sim->p[idx_p] - sim->p[idx_p - sim->ny]);  // dp/dx
            diff = alpha * (RR(u, idx_u, k) + AA(u, idx_u, k) - 4*u[idx_u] + LL(u, idx_u, k) + BB(u, idx_u, k));  // d2u/dx2 + d2u/dy2
            sim->u_star[idx_u] = sim->u[idx_u] + sim->dt * (conv + pres + diff);
        }
    }

    // update v field
    for (i = 1; i <= sim->nx; i++) {
        for (j = 1; j < sim->ny; j++) {
            
            idx_v = i * l + j; 
            idx_p = (i-1) * sim->ny + j;
            
            conv = coef_1 * sim->Hy[idx_v] + coef_2 * sim->Hy_prev[idx_v];  // (u dv/dx + v dv/dy) at t_n and t_(n-1)
            pres = -(sim->p[idx_p] - sim->p[idx_p - 1]);  // dp/dy
            diff = alpha * (RR(v, idx_v, l) + AA(v, idx_v, l) - 4*v[idx_v] + LL(v, idx_v, l) + BB(v, idx_v, l));  // d2v/dx2 + d2v/dy2
            sim->v_star[idx_v] = sim->v[idx_v] + sim->dt * (conv + pres + diff);
        }
    }
}

void update_2(data_Sim *sim) {
    int i, j, idx_u, idx_v, idx_p;
    double dphi;
    int k = sim->ny + 2;
    int l = sim->ny + 1;

    // update u field
    for (i = 1; i < sim->nx; i++) {
        for (j = 1; j <= sim->ny; j++) {
            idx_u = i * k + j;
            idx_p = i * sim->ny + j-1;

            dphi = sim->phi[idx_p] - sim->phi[idx_p-sim->ny];
            sim->u[idx_u] = sim->u_star[idx_u] - sim->dt * dphi;  // - d phi / dx
        }
    }

    // update v field
    for (i = 1; i <= sim->nx; i++) {
        for (j = 1; j < sim->ny; j++) {
            idx_v = i * l + j; 
            idx_p = (i-1) * sim->ny + j;

            dphi = sim->phi[idx_p] - sim->phi[idx_p-1];
            sim->v[idx_v] = sim->v_star[idx_v] - sim->dt * dphi; // - d phi / dy
        }
    }

    // update p field
    for (i = 0; i < sim->size_p; i++) {
        sim->p[i] += sim->phi[i];
    }

    
    // Outflow boundary condition
    /* 
    
    TO BE MODIFIED : REALLY NEED (???) TO KEEP IN MEMORY THESOLUTION AT PREVIOUS TIME ?

    // double w, w_left;
    for (j = 1; j <= sim->ny; j++) {  // u is advected at velocity (1 - uMesh)
        idx_u = sim->nx * k + j;
        sim->u[idx_u] += sim->dt * (1. - sim->uMesh) * (sim->u[idx_u] - LL(sim->u, idx_u, k)) / sim->h;
    }
    for (j = 0; j <= sim->ny; j++) {  // w is advected at velocity (1 - uMesh)
        // get vorticity at the second to last index
        idx_u = (sim->nx - 1) * k + j;  // up left wrt v
        idx_v = (sim->nx + 0) * l + j;
        w_left = (sim->v_prev[idx_v] - LL(sim->v_prev, idx_v, l)) - (sim->u_prev[idx_u] - BB(sim->u_prev, idx_u, k));

        // get vorticity at the last index
        idx_u = (sim->nx + 0) * k + j;  // up left wrt v
        idx_v = (sim->nx + 1) * l + j;
        w = (sim->v_prev[idx_v] - LL(sim->v_prev, idx_v, l)) - (sim->u_prev[idx_u] - BB(sim->u_prev, idx_u, k));

        // derive the condition on v
        sim->v[idx_v] = LL(sim->v, idx_v, l) + (sim->u[idx_u] - sim->u[idx_u-1]) + (1. - sim->uMesh) * sim->dt / sim->h * (w - w_left);
    }*/

}

void integrate_flow(data_Sim *sim, Poisson_data *poisson) {
    int t = 0;

    set_ghost_points(sim);
    set_bd_conditions(sim, sim->u, sim->v);
    set_bd_conditions(sim, sim->u_star, sim->v_star);
    save_fields(sim, t);

    printf("starting ... \n"); fflush(stdout);

    for (t = 1; t <= sim->nt; t++) {

        // sim->uMesh = ALPHA * 2. * M_PI * STROUHAL * sin(2. * M_PI * STROUHAL * t * sim->dt);
        // sim->vMesh = 0.;
        
        compute_convection(sim);
        update_1(sim, t);
        set_bd_conditions(sim, sim->u_star, sim->v_star);
        poisson_solver(sim, poisson);
        update_2(sim);        
        set_bd_conditions(sim, sim->u, sim->v);
        set_ghost_points(sim);

        if (t % 1 == 0) {
            printf("Saving results ... t = %3d\n", t);
            save_fields(sim, t);
        }
    }

    printf("\n");

}


void free_data_Sim(data_Sim *sim) {
    free(sim->u);
    free(sim->v);
    free(sim->p);
}


int main(int argc, char *argv[]){
    // argv : -ksp_type gmres -pc_type lu
    PetscInitialize(&argc, &argv, 0, 0);

    data_Sim *simulation = (data_Sim *)malloc(sizeof(data_Sim));
    init_data_sim(simulation);
    init_fields(simulation);

    Poisson_data *poisson = (Poisson_data *)malloc(sizeof(Poisson_data));
    initialize_poisson_solver(simulation, poisson);

#if DEBUG
    
    poisson_solver(simulation, poisson);
    FILE *ptr = fopen("test_poisson.txt", "w");
    fprintf(ptr, "%d\n", simulation->n);
    for (int i = 0; i < simulation->nx; i++) {
        for (int j = 0; j < simulation->ny; j++) {
            fprintf(ptr, "%.4le ", simulation->phi[i*simulation->ny+j]);
        }
        fprintf(ptr, "\n");
    } 
    fclose(ptr);

#else

    integrate_flow(simulation, poisson);

#endif

    free_data_Sim(simulation);
    free_poisson_solver(poisson);
    PetscFinalize();
}
