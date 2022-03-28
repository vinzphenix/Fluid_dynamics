#include "poisson.h"
#include "project.h"

char *myPath = "./data/";
char filename_params[50];
char filename_u[50];
char filename_v[50];
char filename_p[50];


void init_data_sim(data_Sim *sim) {
    sprintf(filename_params, "%ssimu_params.txt", myPath);
    sprintf(filename_u, "%ssimu_u.txt", myPath);
    sprintf(filename_v, "%ssimu_v.txt", myPath);
    sprintf(filename_p, "%ssimu_p.txt", myPath);

    // sim->h = CFL / (FOURIER * RE * 5);
    // sim->dt = FOURIER * RE * (sim->h) * (sim->h);

    // sim->nt = ceil(TEND / sim->dt);
    // sim->nx = (double) L_ / sim->h;
    // sim->ny = (double) H_ / sim->h;

    sim->n = 1;  // 1. / h

    sim->nt = 1;
    sim->nx = 15 * sim->n;
    sim->ny = 5 * sim->n;

    sim->dt = TEND / (double) sim->nt;
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
    sim->u_temp = sim->u + size_u;
    sim->u_prev = sim->u + 2*size_u;
    sim->Hx = sim->u + 3*size_u;
    sim->Hx_prev = sim->u + 4*size_u;

    sim->v = (double *)calloc(5*size_v, sizeof(double));
    sim->v_temp = sim->v + size_v;
    sim->v_prev = sim->v + 2*size_v;
    sim->Hy = sim->v + 3*size_v;
    sim->Hy_prev = sim->v + 4*size_v;

    sim->p = (double *)calloc(2*size_p, sizeof(double));
    sim->phi = sim->p + size_p;
}


void save_fields(data_Sim *sim, int t) {
    FILE *ptr, *ptr_u, *ptr_v, *ptr_p;
    if (t == 0) {
        ptr = fopen(filename_params, "w");
        fprintf(ptr, "%d %d %d\n", sim->nt, sim->nx, sim->ny);
        fprintf(ptr, "%lf %lf %lf %lf %lf %lf %lf %lf\n", TEND, L_, H_, LBOX, D_IN, D_BOT, sim->dt, sim->h);
    } else {
        ptr_u = fopen(filename_u, "a");
        ptr_v = fopen(filename_v, "a");
        ptr_p = fopen(filename_p, "a");
    }

    for (int i = 0; i < sim->size_u; i++) {
        fprintf(ptr, "%.5le ", sim->u[i]);
    }
    fprintf(ptr, "\n");

    for (int i = 0; i < sim->size_v; i++) {
        fprintf(ptr, "%.5le ", sim->v[i]);
    }
    fprintf(ptr, "\n");

    for (int i = 0; i < sim->size_p; i++) {
        fprintf(ptr, "%.5le ", sim->p[i]);
    }
    fprintf(ptr, "\n");
    
    if (t == 0) {
        fclose(ptr);
    } else {
        fclose(ptr_u);
        fclose(ptr_v);
        fclose(ptr_p);
    }
}


void compute_shadow_pts(data_Sim *sim) {
    int i, j;
    // int idx_u, idx_v;
    
    int k = sim->ny + 2;
    int l = sim->ny + 1;

    // Inflow boundary condition
    for (j = 1; j < sim->ny; j++) {   // Inflow uniform profile u
        sim->u[j] = 1.;               // could do a special profile
    }
    for (j = 0; j <= sim->ny; j++) {  // trick for zero v at inflow
        sim->v[j] = -0.2 * (sim->v[3*l + j] - 5. * sim->v[2*l + j] + 15. * sim->v[l + j] - 16. * 0.);
    }

    // external walls (v = 0 AND w = 0 --> du/dy = 0)
    for (i = 0; i < sim->nx + 1; i++) {                // zero voriticty
        sim->u[i*k] = sim->u[i*k + 1];                 // below
        sim->u[i*k + (k-1)] = sim->u[i*k + (k-2)];     // above
    }
    for (i = 0; i < sim->nx + 2; i++) { // no-through flow
        sim->v[i*l] = 0.;               // below
        sim->v[i*l + (l-1)] = 0.;       // above
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


void init_fields(data_Sim *sim) {
    int i;
    srand(12102000);
    for (i = 0; i < sim->size_u; i++) {
        sim->u[i] = 0 + (rand() % (20 - 0 + 1));
    }
    for (i = 0; i < sim->size_v; i++) {
        sim->v[i] = 0 + (rand() % (20 - 0 + 1));
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


void update_1(data_Sim *sim) {
    int i, j, idx_u, idx_v, idx_p;
    double conv, pres, diff;
    
    int k = sim->ny + 2;
    int l = sim->ny + 1;

    double *u = sim->u;
    double *v = sim->v;

    // update u field
    for (i = 1; i < sim->nx; i++) {
        for (j = 1; j <= sim->ny; j++) {
            
            idx_u = i * k + j;
            idx_p = i * sim->ny + j-1;
            
            conv = 0.5 * (3. * sim->Hx[idx_u] - sim->Hx_prev[idx_u]);  // u du/dx + v du/dy
            pres = sim->p[idx_p] - sim->p[idx_p - sim->ny];  // dp/dx
            diff = 1. / (RE * pow(sim->h, 2)) * (RR(u, idx_u, k) + AA(u, idx_u, k) - 4*u[idx_u] + LL(u, idx_u, k) + BB(u, idx_u, k));  // d2u/dx2 + d2u/dy2
            sim->u_temp[idx_u] = sim->u[idx_u] + sim->dt * (-conv - pres + diff);
        }
    }

    // update v field
    for (i = 1; i <= sim->nx; i++) {
        for (j = 1; j < sim->ny; j++) {
            
            idx_v = i * l + j; 
            idx_p = (i-1) * sim->ny + j;
            
            conv = 0.5 * (3. * sim->Hy[idx_v] - sim->Hy_prev[idx_v]);  // u dv/dx + v dv/dy
            pres = sim->p[idx_p] - sim->p[idx_p - 1];  // dp/dy
            diff = 1. / (RE * pow(sim->h, 2)) * (RR(v, idx_v, l) + AA(v, idx_v, l) - 4*v[idx_v] + LL(v, idx_v, l) + BB(v, idx_v, l));  // d2v/dx2 + d2v/dy2
            sim->v_temp[idx_v] = sim->v[idx_v] + sim->dt * (-conv - pres + diff);
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
            sim->u[idx_u] = sim->u_temp[idx_u] - sim->dt * dphi;  // - d phi / dx
        }
    }

    // update v field
    for (i = 1; i <= sim->nx; i++) {
        for (j = 1; j < sim->ny; j++) {
            idx_v = i * l + j; 
            idx_p = (i-1) * sim->ny + j;

            dphi = sim->phi[idx_p] - sim->phi[idx_p-1];
            sim->v[idx_v] = sim->v_temp[idx_v] - sim->dt * dphi; // - d phi / dy
        }
    }

    // update p field
    for (i = 0; i < sim->size_p; i++) {
        sim->p[i] += sim->phi[i];
    }

    double w, w_left;
    // Outflow boundary condition
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
    }

}

void integrate_flow(data_Sim *sim) {
    for (int t = 1; t <= sim->nt; t++) {
        
        //compute_shadow_pts(sim);
        // compute_convection(sim);
        // update_1(sim);
        // poisson ???
        // update_2(sim);

    }

    // compute_convection(sim);
    // save_debug(sim);

}

int main(int argc, char *argv[]){
    
    PetscInitialize(&argc, &argv, 0, 0);

    data_Sim *simulation = (data_Sim *)malloc(sizeof(data_Sim));
    init_data_sim(simulation);
    init_fields(simulation);

    integrate_flow(simulation);

    PetscFinalize();

}
