#include "poisson.h"
#include "project.h"
#include "adi.h"

char filename_params[50];
char filename_stats[50];
char filename_u[50];
char filename_v[50];
char filename_p[50];
char filename_T[50];
char filename_u_avg[50];
char filename_v_avg[50];

void init_Sim_data(Sim_data *sim, int n_input, double dt_input, double tend_input, double freq, const char *myPath) {

    sprintf(filename_params, "%ssimu_params.txt", myPath);
    sprintf(filename_stats, "%ssimu_stats.txt", myPath);
    sprintf(filename_u, "%ssimu_u.txt", myPath);
    sprintf(filename_v, "%ssimu_v.txt", myPath);
    sprintf(filename_p, "%ssimu_p.txt", myPath);
    sprintf(filename_T, "%ssimu_T.txt", myPath);
    sprintf(filename_u_avg, "%ssimu_u_avg.txt", myPath);
    sprintf(filename_v_avg, "%ssimu_v_avg.txt", myPath);
    
    // erase from previous run
    fclose(fopen(filename_stats, "w"));
    fclose(fopen(filename_u, "w"));
    fclose(fopen(filename_v, "w"));
    fclose(fopen(filename_p, "w"));
    fclose(fopen(filename_T, "w"));

    // Space discretization
    sim->n = n_input;
    sim->nx = L_ * sim->n;
    sim->ny = H_ * sim->n;
    sim->h = 1. / ((double) sim->n);

    // Time discretization
    sim->dt_stable = fmin(FOURIER * RE * (sim->h) * (sim->h), CFL * sim->h / U0V0);
    sim->dt = (dt_input == 0.) ? sim->dt_stable : dt_input;
    sim->tsim = tend_input;
    sim->nt = (int) ceil(sim->tsim / sim->dt);
    printf("Current \u0394t = %.5lf  vs  %.5lf = \u0394t maximal\n", sim->dt, sim->dt_stable);

    sim->tnow = 0.;
    sim->elapsed = 0.;
    sim->save_freq = freq;
    if (START_AVG > sim->tsim - 1) {
        sim->start_avg = 0.75 * sim->tsim;
    } else {
        sim->start_avg = START_AVG;
    }
    
    // Obstacle initially not moving
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
    

    // Allocate u and set matrix access
    sim->u_data = (double *)calloc(5 * size_u, sizeof(double));
    sim->u_avg = sim->u_data + 4 * size_u;

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
    sim->v_data = (double *)calloc(5 * size_v, sizeof(double));
    sim->v_avg = sim->v_data + 4 * size_v;
    
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


#   if TEMP_MODE
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

#   if TEMP_MODE
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


    // set u to U_inf at the beginning (0 is already good for v)
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
    for (j = 1; j < sim->ny + 1; j++) {
        eta = (j-0.5) * sim->h / H_;
        sim->U[0][j] = 6. * 1. * eta * (1. - eta);     // Poiseuille
        sim->US[0][j] = sim->U[0][j];                  // Poiseuille
    }
    for (i = 0; i < sim->nx + 1; i++) {
        j = 0;
        sim->U[i][j] = -0.2 * (sim->U[i][j+3] - 5.*sim->U[i][j+2] + 15.*sim->U[i][j+1] - 16 * 0.);
        j = sim->ny+1;
        sim->U[i][j] = -0.2 * (sim->U[i][j-3] - 5.*sim->U[i][j-2] + 15.*sim->U[i][j-1] - 16 * 0.);
    }
#   endif

#   if WALL_DIRICHLET
    double one_to_zero;
    for (i = 1; i < sim->nx + 1; i++) {
        // one_to_zero =  0.5 * (1. - tanh( ((i - 0.5) * sim->h - 0.75 * L_) / 1.));
        one_to_zero = 1.;
        j = 0;
        sim->T[i][j] = -0.2 * (sim->T[i][j+3] - 5.*sim->T[i][j+2] + 15.*sim->T[i][j+1] - 16.*TMAX*one_to_zero);
        j = sim->ny+1;
        sim->T[i][j] = -0.2 * (sim->T[i][j-3] - 5.*sim->T[i][j-2] + 15.*sim->T[i][j-1] - 16.*TMIN*one_to_zero);
    }
#   endif

#   if BOX_BOT_TOP_DIRICHLET
    for (i = sim->n * D_IN + 1; i < sim->n * (D_IN + LBOX) + 1; i++) {
        j = sim->n * D_BOT + 1;  // lower wall of rectangle
        sim->T[i][j] = -0.2 * (sim->T[i][j-3] - 5.*sim->T[i][j-2] + 15.*sim->T[i][j-1] - 16.*TMAX);
        j = sim->n * (D_BOT + 1);  // upper wall of rectangle
        sim->T[i][j] = -0.2 * (sim->T[i][j+3] - 5.*sim->T[i][j+2] + 15.*sim->T[i][j+1] - 16.*TMIN);
    }
#   endif

#   if BOX_LFT_RGT_DIRICHLET
    i = sim->n * D_IN + 1;  // left wall of rectanle
    for (j = sim->n * D_BOT + 1; j < sim->n * (D_BOT + 1) + 1; j++)
        sim->T[i][j] = -0.2 * (sim->T[i-3][j] - 5.*sim->T[i-2][j] + 15.*sim->T[i-1][j] - 16.*TMAX);
    i = sim->n * (D_IN + LBOX);  // right wall of rectanle
    for (j = sim->n * D_BOT + 1; j < sim->n * (D_BOT + 1) + 1; j++)
        sim->T[i][j] = -0.2 * (sim->T[i+3][j] - 5.*sim->T[i+2][j] + 15.*sim->T[i+1][j] - 16.*TMAX);
#   endif
}


int increase_dt = 0; 
int decrease_dt = 0;
// int nChanges = 0;
double dt_change;
double dt_init;
void adapt_dt(Sim_data *sim, double speed) {

    // Handle adaptative timestep
    if (sim->tnow < 1.) {
        dt_init = sim->dt;
        return;
    }
    int duration = 10;

    if ((decrease_dt == 0) && (increase_dt == 0)){
        // keep dt in [0.90, 0.95]*dt_stable
        if ((sim->dt > sim->dt_stable * 0.95) && (sim->dt > dt_init * 0.1)) {
            // second condition avoids doing smaller and smaller time steps
            // ((sim->rew > 40.) || (sim->reh > 40.)) ||
            decrease_dt = 1;
            dt_change = sim->dt_stable * 0.05;
            // nChanges --;  // not used anymore
        } else if ((sim->dt < sim->dt_stable * 0.90))  {
            //((sim->rew < 32.) && (sim->reh < 32.)) || 
            increase_dt = 1;
            dt_change = sim->dt_stable * 0.05;
            // nChanges ++;  // not used anymore
        }
    }
    if ((0 < decrease_dt) && (decrease_dt <= duration)) {
        sim->dt -= dt_change / duration;
        decrease_dt ++;
    } else if ((0 < increase_dt) && (increase_dt <= duration)) {
        sim->dt += dt_change / duration;
        increase_dt ++;
    } else if (decrease_dt > duration) {
        decrease_dt = 0;
    } else if (increase_dt > duration) {
        increase_dt = 0;
    }
}


void save_fields(Sim_data *sim) {
    FILE *ptr, *ptr_u, *ptr_v, *ptr_p, *ptr_T;
    int i;

    if (sim->tnow == 0) {
        ptr = fopen(filename_params, "w");
        fprintf(ptr, "%d %d %d %d\n", sim->nx, sim->ny, sim->n, TEMP_MODE);
        fprintf(ptr, "%le %le %.10le %le %le %le\n", RE, sim->tsim, sim->dt, sim->h, sim->save_freq, sim->start_avg);
        fprintf(ptr, "%d %d %d %d %d\n", L_, H_, LBOX, D_IN, D_BOT);
        fprintf(ptr, "%le %le %le %le %le %le %le %d\n", ALPHA, STROUHAL, SIWNG_START, KAPPA_Y, STROUHAL_Y, PERT_START, SMOOTH, N_CYCLES);
        fprintf(ptr, "%le %le %le %le %le\n", PR, GR, EC, TMIN, TMAX);
        fclose(ptr);
    }

    ptr_u = fopen(filename_u, "a");
    ptr_v = fopen(filename_v, "a");
    ptr_p = fopen(filename_p, "a");
    ptr_T = fopen(filename_T, "a");

    for (i = 0; i < sim->size_u; i++) fprintf(ptr_u, FMT, sim->u_data[i]);  // velocity x
    for (i = 0; i < sim->size_v; i++) fprintf(ptr_v, FMT, sim->v_data[i]);  // velocity y
    for (i = 0; i < sim->size_p; i++) fprintf(ptr_p, FMT, sim->p_data[i]);  // pressure

#   if TEMP_MODE
    int j;
    for (i = 1; i < sim->nx + 1; i++) {
        for (j = 1; j < sim->ny + 1; j++) {
            fprintf(ptr_T, FMT, sim->T[i][j]);  // temperature (do not save ghosts)
        }
    }
#   endif

    //close
    fclose(ptr_u);
    fclose(ptr_v);
    fclose(ptr_p);
    fclose(ptr_T);
 }


void save_diagnostics(Sim_data *sim, int saving) {
    
    int i, j, i1, i2, j1, j2;
    double uij, vij, wij;
    double reh = 0.;
    double rew = 0.;
    double drag_p = 0.;
    double lift_p = 0.;
    double drag_f = 0.;
    double lift_f = 0.;
    FILE *ptr;

    // Update averaged fields ("as soon as vortex shedding established")
    if (sim->tnow > sim->start_avg) {
        for (i = 0; i < sim->size_u; i++) sim->u_avg[i] += (sim->u_data[i]) * sim->dt;
        for (i = 0; i < sim->size_v; i++) sim->v_avg[i] += (sim->v_data[i]) * sim->dt;
    }
    
    if (sim->tnow + 1.e-9 >= sim->tsim) {  // save at last iteration
        ptr = fopen(filename_u_avg, "w");
        for (i = 0; i < sim->size_u; i++) fprintf(ptr, "%.10le\n", sim->u_avg[i] / (sim->tnow - sim->start_avg));
        fclose(ptr);

        ptr = fopen(filename_v_avg, "w");
        for (i = 0; i < sim->size_v; i++) fprintf(ptr, "%.10le\n", sim->v_avg[i] / (sim->tnow - sim->start_avg));
        fclose(ptr);
    }


    // Mesh reynolds numbers
    for (i = 0; i < sim->nx+1; i++) {
        for (j = 0; j < sim->ny+1; j++) {
            uij = 0.5 * (sim->U[i][j+1] + sim->U[i][j]);
            vij = 0.5 * (sim->V[i+1][j] + sim->V[i][j]);
            reh = fmax(reh, fabs(uij) + fabs(vij));

            wij = (sim->V[i+1][j] - sim->V[i][j]) - (sim->U[i][j+1] - sim->U[i][j]);
            rew = fmax(rew, fabs(wij / sim->h));
        }
    }
    sim->dt_stable = fmin(FOURIER * RE * (sim->h) * (sim->h), CFL * sim->h / reh);

    // Aerodynamic forces
    i1 = D_IN * sim->n - 1;   i2 = (D_IN + LBOX) * sim->n;
    j1 = D_BOT * sim->n - 1;  j2 = (D_BOT + 1) * sim->n;

    j = j1;
    drag_p += 0.5 * (sim->P[i1][j] - sim->P[i2][j]);
    for (j++; j < j2; j++)
        drag_p += sim->P[i1][j] - sim->P[i2][j];
    drag_p += 0.5 * (sim->P[i1][j] - sim->P[i2][j]);

    i = i1;
    lift_p += 0.5 * (sim->P[i][j1] - sim->P[i][j2]);
    for (i++; i < i2; i++)
        lift_p += sim->P[i][j1] - sim->P[i][j2];
    lift_p += 0.5 * (sim->P[i][j1] - sim->P[i][j2]);


    i1 = D_IN * sim->n;   i2 = (D_IN + LBOX) * sim->n + 1;
    j1 = D_BOT * sim->n;  j2 = (D_BOT + 1) * sim->n + 1;

    j = j1;
    lift_f += 0.5 * ((sim->V[i1][j] - sim->V[i1+1][j]) + (sim->V[i2][j] - sim->V[i2-1][j]));  // trapezoidal integration : first
    for (j++; j < j2 - 1; j++)
        lift_f += (sim->V[i1][j] - sim->V[i1+1][j]) + (sim->V[i2][j] - sim->V[i2-1][j]);      // [ -(dv/dx) left + (dv/dx) right ] * dy
    lift_f += 0.5 * ((sim->V[i1][j] - sim->V[i1+1][j]) + (sim->V[i2][j] - sim->V[i2-1][j]));  // trapezoidal integration : last

    i = i1;
    drag_f += 0.5 * ((sim->U[i][j1] - sim->U[i][j1+1]) + (sim->U[i][j2] - sim->U[i][j2-1]));  // trapezoidal integration : first
    for (i++; i < i2 - 1; i++)
        drag_f += (sim->U[i][j1] - sim->U[i][j1+1]) + (sim->U[i][j2] - sim->U[i][j2-1]);      // [ -(du/dy) below + (du/dy) above ] * dx
    drag_f += 0.5 * ((sim->U[i][j1] - sim->U[i][j1+1]) + (sim->U[i][j2] - sim->U[i][j2-1]));  // trapezoidal integration : last
    
    sim->reh = reh * sim->h * RE;
    sim->rew = rew * sim->h * sim->h * RE;
    drag_p = 2 * drag_p * sim->h;
    drag_f = 2 * drag_f / RE;
    lift_p = 2 * lift_p * sim->h;
    lift_f = 2 * lift_f / RE;
    // factor 2 because Cd = Fd / (0.5 * rho * u^2 * Hbox)

    ptr = fopen(filename_stats, "a");
    fprintf(ptr, "%le %le %le %le %le %le %.10le %d\n", sim->reh, sim->rew, drag_p, drag_f, lift_p, lift_f, sim->tnow, saving);
    fclose(ptr);

#   if ADAPTATIVE_DT
    adapt_dt(sim, reh);
#   endif
}


void set_bd_conditions(Sim_data *sim) {
    int i, j;
    double coef;
    
    /*// Inflow boundary condition for u  // never changes
    for (j = 1; j <= sim->ny; j++) {
        U[0][j] = 1.;
    }*/

    // Outflow boundary condition for u
    i = sim->nx;
    coef = sim->dt / sim->h * (1. - sim->uMesh);
    for (j = 1; j <= sim->ny; j++) {  // u is advected at velocity (1 - uMesh)
        sim->US[i][j] = sim->U[i][j] - coef * (sim->U[i][j] - sim->U[i-1][j]);
        sim->U[i][j] = sim->US[i][j];
    }

#   if NO_SLIP  // if slip: v=0 at all times
    // Channel walls: condition for v
    for (i = 0; i < sim->nx + 2; i++) { // no-through flow
        sim->VS[i][0] = sim->vMesh;             // below
        sim->VS[i][sim->ny] = sim->vMesh;       // above
        
        sim->V[i][0] = sim->VS[i][0];
        sim->V[i][sim->ny] = sim->VS[i][sim->ny];
    }
#   endif


    i = D_IN * sim->n;  // Left wall of rectangle (u)
    for (j = D_BOT * sim->n + 1; j <= (D_BOT + 1) * sim->n; j++) {
        sim->US[i][j] = sim->uMesh;
        sim->U[i][j] = sim->US[i][j];
    }
    
    i = (D_IN + LBOX) * sim->n; // Right wall of rectangle (u)
    for (j = D_BOT * sim->n + 1; j <= (D_BOT + 1) * sim->n; j++) {
        sim->US[i][j] = sim->uMesh;
        sim->U[i][j] = sim->US[i][j];
    }


    for (i = D_IN * sim->n + 1; i < (D_IN + LBOX) * sim->n; i++) {                
        
        j = D_BOT * sim->n;  // lower wall of the rectangle
        sim->VS[i][j] = sim->vMesh;
        sim->V[i][j] = sim->VS[i][j];

        j = (D_BOT + 1) * sim->n;  // upper wall of the rectangle
        sim->VS[i][j] = sim->vMesh;
        sim->V[i][j] = sim->VS[i][j];
    }
}

void set_ghost_points(Sim_data *sim) {
    int i, j;
    double coef, w_last, w_left;

    // External walls
    for (i = 0; i < sim->nx + 1; i++) {
        j = 0;
#       if NO_SLIP  // (v = 0) and (u = 0)
            sim->U[i][j] = -0.2 * (sim->U[i][j+3] - 5.*sim->U[i][j+2] + 15.*sim->U[i][j+1] - 16.*sim->uMesh);
#       else
            sim->U[i][0] = sim->U[i][1];  // (v = 0 and w = 0 --> du/dy = 0)
#       endif

        j = sim->ny+1;
#       if NO_SLIP
            sim->U[i][j] = -0.2 * (sim->U[i][j-3] - 5.*sim->U[i][j-2] + 15.*sim->U[i][j-1] - 16.*sim->uMesh);
#       else
            sim->U[i][j] = sim->U[i][j-1];
#       endif
    }


    // Inflow boundary condition for v
    for (j = 1; j < sim->ny; j++) {  // trick with ghost for zero v at inflow
        sim->V[0][j] = -0.2 * (sim->V[3][j] - 5. * sim->V[2][j] + 15. * sim->V[1][j] - 16. * 0.);  // 0 or vMesh ?
    }
    
    // Outflow boundary condition for v
    // write into VS, and then copy back into V to avoid using updated values V(n+1) in the computation of V(n+1)
    i = sim->nx;  // before last i_index of v_ij
    coef = sim->dt * (1. - sim->uMesh);
    for (j = 1; j < sim->ny; j++) {
        w_last = ((sim->V[i+1][j] - sim->V[i  ][j]) - (sim->U[i  ][j+1] - sim->U[i  ][j])) / sim->h;
        w_left = ((sim->V[i  ][j] - sim->V[i-1][j]) - (sim->U[i-1][j+1] - sim->U[i-1][j])) / sim->h;
        sim->VS[i+1][j] = (sim->V[i][j] + sim->U[i][j+1] - sim->U[i][j]) + sim->h * w_last - coef * (w_last - w_left);
        // [(sim->VS[i+1][j] - sim->V[i][j]) - (sim->U[i][j+1] - sim->U[i][j])/h - w_last ]/dt = -uc * (w_last - w_left) / dx;
    }
    for (j = 1; j < sim->ny; j++) {
        sim->V[i+1][j] = sim->VS[i+1][j];
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
    double coef;

    for (i = 1; i < sim->nx + 1; i++) {
#       if WALL_DIRICHLET
            // coef =  0.5 * (1. - tanh( ((i - 0.5) * sim->h - 0.75 * L_) / 1.));
            coef = 1.;
            j = 0;  // Below
            sim->T[i][j] = -0.2 * (sim->T[i][j+3] - 5.*sim->T[i][j+2] + 15.*sim->T[i][j+1] - 16. * TMAX * coef);  // Dirichlet
            j = sim->ny+1;  // Above
            sim->T[i][j] = -0.2 * (sim->T[i][j-3] - 5.*sim->T[i][j-2] + 15.*sim->T[i][j-1] - 16. * TMIN * coef);  // Dirichlet
#       else
            j = 0;  // Below
            sim->T[i][j] = sim->T[i][j+1];  // Adiabatic
            j = sim->ny+1;  // Above
            sim->T[i][j] = sim->T[i][j-1];  // Adiabatic
#       endif
    }

    coef = sim->dt / sim->h * (1. - sim->uMesh);
    for (j = 1; j < sim->ny + 1; j++) sim->T[0        ][j] = sim->T[1][j];                                         // inflow
    for (j = 1; j < sim->ny + 1; j++) sim->T[sim->nx+1][j] -= coef * (sim->T[sim->nx+1][j] - sim->T[sim->nx][j]);  // outflow (advected ?)

    i = sim->n * D_IN + 1;  // left wall of rectanle
    for (j = sim->n * D_BOT + 1; j < sim->n * (D_BOT + 1) + 1; j++) {
#       if BOX_LFT_RGT_DIRICHLET
        sim->T[i][j] = -0.2 * (sim->T[i-3][j] - 5.*sim->T[i-2][j] + 15.*sim->T[i-1][j] - 16.*TMAX);  // Dirichlet
#       else
        sim->T[i][j] = sim->T[i-1][j];  // No flux
#       endif
    }

    i = sim->n * (D_IN + LBOX);  // right wall of rectanle
    for (j = sim->n * D_BOT + 1; j < sim->n * (D_BOT + 1) + 1; j++) {
#       if BOX_LFT_RGT_DIRICHLET
        sim->T[i][j] = -0.2 * (sim->T[i+3][j] - 5.*sim->T[i+2][j] + 15.*sim->T[i+1][j] - 16.*TMAX);  // Dirichlet
#       else
        sim->T[i][j] = sim->T[i+1][j];  // No flux
#       endif
    }
    // if dirichlet applied on lateral walls of the box --> change special treatment in corrector step
    // sim->T[i][j] = -0.2 * (sim->T[i+3][j] - 5.*sim->T[i+2][j] + 15.*sim->T[i+1][j] - 16.*TMAX);


    // ! Ghost points in the corners are overwritten
    for (i = sim->n * D_IN + 1; i < sim->n * (D_IN + LBOX) + 1; i++) {
        
        j = sim->n * D_BOT + 1;  // lower wall of rectangle
#       if BOX_BOT_TOP_DIRICHLET
            sim->T[i][j] = -0.2 * (sim->T[i][j-3] - 5.*sim->T[i][j-2] + 15.*sim->T[i][j-1] - 16.*TMAX);
#       else
            sim->T[i][j] = sim->T[i][j-1];
#       endif

        j = sim->n * (D_BOT + 1);  // upper wall of rectangle
#       if BOX_BOT_TOP_DIRICHLET
            sim->T[i][j] = -0.2 * (sim->T[i][j+3] - 5.*sim->T[i][j+2] + 15.*sim->T[i][j+1] - 16.*TMIN);
#       else
            sim->T[i][j] = sim->T[i][j+1];
#       endif
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

    double **U = sim->U;
    double **V = sim->V;

    double coef_1 = (sim->first_iteration) ? -1. : -1.5;
    double coef_2 = (sim->first_iteration) ? +0. : +0.5;
    double nu_h2 = 1. / (RE * sim->h * sim->h);
#   if TEMP_MODE
    double beta = 0.5 * GR / (RE * RE);  // 0.5 here to prepare the average in the "for" loop
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

                sim->US[i][j] = U[i][j] + sim->dt * (
                    + coef_1 * sim->HX[i][j] + coef_2 * sim->HX_[i][j]                     // Convection
                    - (sim->P[i][j-1] - sim->P[i-1][j-1]) * sim->n                         // Pressure gradient
                    + nu_h2 * (U[i+1][j] + U[i-1][j] - 4*U[i][j] + U[i][j+1] + U[i][j-1])  // Diffusion
                );
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

                sim->VS[i][j] = V[i][j] + sim->dt * (
                    + coef_1 * sim->HY[i][j] + coef_2 * sim->HY_[i][j]                    // Convection
                    - (sim->P[i-1][j] - sim->P[i-1][j-1]) * sim->n                        // Pressure gradient
                    + nu_h2 * (V[i+1][j] + V[i-1][j] - 4*V[i][j] + V[i][j+1] + V[i][j-1]) // Diffusion
#                   if TEMP_MODE
                    + beta * (sim->T[i][j+1] + sim->T[i][j])  // AVG(T above, T below)  // + sign because g has negative y component
#                   endif
                );
            }
        }
    }
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
}


void handle4corners(Sim_data *sim, double diff) {
    // Handle 4 points where ghost points were overwritten (corners of the rectangle)
    int i, j;
    
#   if BOX_LFT_RGT_DIRICHLET
    double ghost;
    i = (D_IN) * sim->n;
    j = (D_BOT) * sim->n + 1;
    ghost = -0.2 * (sim->T[i-2][j] - 5.*sim->T[i-1][j] + 15.*sim->T[i][j] - 16.*TMAX);  // Set ghost
    sim->PHI[i-1][j-1] += diff * (ghost - sim->T[i+1][j]);  // Add good ghost and delete bad one
    j = (D_BOT + 1) * sim->n;
    ghost = -0.2 * (sim->T[i-2][j] - 5.*sim->T[i-1][j] + 15.*sim->T[i][j] - 16.*TMAX);
    sim->PHI[i-1][j-1] += diff * (ghost - sim->T[i+1][j]);
    
    i = (D_IN + LBOX) * sim->n + 1;
    j = (D_BOT) * sim->n + 1;
    ghost = -0.2 * (sim->T[i+2][j] - 5.*sim->T[i+1][j] + 15.*sim->T[i][j] - 16.*TMAX);
    sim->PHI[i-1][j-1] += diff * (ghost - sim->T[i-1][j]);
    j = (D_BOT + 1) * sim->n;
    ghost = -0.2 * (sim->T[i+2][j] - 5.*sim->T[i+1][j] + 15.*sim->T[i][j] - 16.*TMAX);
    sim->PHI[i-1][j-1] += diff * (ghost - sim->T[i-1][j]);

#   else
    i = (D_IN) * sim->n;
    j = (D_BOT) * sim->n + 1;
    sim->PHI[i-1][j-1] += diff * (sim->T[i][j] - sim->T[i+1][j]);  // T[i+1][j] was the ghost for the temperature below the lower-left corner of the box
    j = (D_BOT + 1) * sim->n;
    sim->PHI[i-1][j-1] += diff * (sim->T[i][j] - sim->T[i+1][j]);  // T[i+1][j] was the ghost for the temperature above the upper-left corner of the box
    i = (D_IN + LBOX) * sim->n + 1;
    j = (D_BOT) * sim->n + 1;
    sim->PHI[i-1][j-1] += diff * (sim->T[i][j] - sim->T[i-1][j]);  // T[i-1][j] was the ghost for the temperature below the lower-right corner of the box
    j = (D_BOT + 1) * sim->n;
    sim->PHI[i-1][j-1] += diff * (sim->T[i][j] - sim->T[i-1][j]);  // T[i-1][j] was the ghost for the temperature above the upper-right corner of the box

#   endif

}

void corrector_step_temperature(Sim_data *sim) {
    int i, j;
    int i_s, i_f, j_s, j_f;

    double coef_1 = (sim->first_iteration) ? -1. : -1.5;
    double coef_2 = (sim->first_iteration) ? +0. : +0.5;
    double diff = 1. / (PR * RE) * 1. / (sim->h * sim->h);
    
    double diss = (EC / RE) * 1. / (sim->h * sim->h);
    double dudx, cross, dvdy;

    // for (i = 1; i < sim->nx + 1; i++) {
        // for (j = 1; j < sim->ny + 1; j++) {

    // compute rhs of temperature equation (diffusion + viscous dissipation)
    for (int block = 8; block < 12; block++) {
        i_s = sim->i_start[block];
        i_f = sim->i_final[block];
        j_s = sim->j_start[block];
        j_f = sim->j_final[block];

        for (i = i_s; i < i_f; i++) {
            for (j = j_s; j < j_f; j++) {
                dudx = (sim->U[i][j] - sim->U[i-1][j]);  // denominator is dx (already in diss)
                dvdy = (sim->V[i][j] - sim->V[i][j-1]);  // denominator is dy (already in diss)
                cross = 0.25 * ((sim->U[i-1][j+1] + sim->U[i][j+1]) - (sim->U[i-1][j-1] + sim->U[i][j-1])    // denominator is 2*2dx
                            + (sim->V[i+1][j-1] + sim->V[i+1][j]) - (sim->V[i-1][j-1] + sim->V[i-1][j]));  // denominator is 2*2dy
                
                sim->PHI[i-1][j-1] =     // use existing array PHI for this purpose -> index -1 since size of PHI smaller
                    diff * (sim->T[i+1][j] + sim->T[i-1][j] + sim->T[i][j+1] + sim->T[i][j-1] - 4. * sim->T[i][j])
                  + diss * (2. * (dudx * dudx + dvdy * dvdy) + cross * cross);
            }
        }
    }

    // ! The problem is that the 4 corners inside the rectangle are ghost points of 2 different cells each
    handle4corners(sim, diff);  // Does it really make a difference ? Not sure, since ghost values should not be very much different

    for (int block = 8; block < 12; block++) {
        i_s = sim->i_start[block];
        i_f = sim->i_final[block];
        j_s = sim->j_start[block];
        j_f = sim->j_final[block];

        for (i = i_s; i < i_f; i++) {
            for (j = j_s; j < j_f; j++) {
                sim->T[i][j] += sim->dt * (coef_1 * sim->HT[i][j] + coef_2 * sim->HT_[i][j] + sim->PHI[i-1][j-1]);
                // PHI is just an array that contains both DIFFUSION and DISSIPATION
            }
        }
    }
}


void swap_next_previous(Sim_data *sim) {
    double **tmp;
    
    tmp = sim->HX;
    sim->HX = sim->HX_;
    sim->HX_ = tmp;

    tmp = sim->HY;
    sim->HY = sim->HY_;
    sim->HY_ = tmp;

#   if TEMP_MODE
    tmp = sim->HT;
    sim->HT = sim->HT_;
    sim->HT_ = tmp;
#   endif
}


void set_mesh_velocity(Sim_data *sim, double t_now) {
    // double t_now = sim->dt * sim->t;
    double pert_end = PERT_START + ((double) N_CYCLES / STROUHAL_Y);
    
    if (SIWNG_START < t_now) {
        sim->uMesh = ALPHA * sin(2. * M_PI * STROUHAL * (t_now - SIWNG_START));
    } else {
        sim->uMesh = 0.;
    }

    if ((PERT_START <= t_now) && (t_now <= PERT_START + ((double) N_CYCLES / STROUHAL_Y))) {
        sim->vMesh = KAPPA_Y * 2. * M_PI * STROUHAL_Y * sin(2. * M_PI * STROUHAL_Y * (t_now - PERT_START));

        if (SMOOTH > 1e-3) {  // otherwise, no smoothing
            double f, g1, g2, dg1, dg2;  // f is position, g1 and g2 are smoothing functions at start and end, dg1 and dg2 are derivatives
            f = KAPPA_Y * (1. - cos(2. * M_PI * STROUHAL_Y * (t_now - PERT_START)));
            g1 = 1. - exp((PERT_START - t_now) / SMOOTH);
            g2 = 1. - exp((t_now - pert_end) / SMOOTH);
            dg1 = 1 / SMOOTH * exp((PERT_START - t_now) / SMOOTH);
            dg2 = -1. / SMOOTH * exp((t_now - pert_end) / SMOOTH);
            sim->vMesh = sim->vMesh * g1 * g2 + f * (dg1 * g2 + g1 * dg2);
        }

    } else {
        sim->vMesh = 0.;
    }
}


void check_boundary(Sim_data *sim) {
    int i, j;

    // set identity for upper and lower walls
    for (i = 0; i < sim->nx + 1; i++) {
        for (j = 0; j < sim->ny + 2; j += sim->ny+1) {
            printf("u[%d][%d] = %lf \t us = %lf --> Upper an lower\n", i, j, sim->U[i][j], sim->US[i][j]);
        }
    }

    // set identity for left and right walls
    for (j = 1; j < sim->ny + 1; j++) {
        for (i = 0; i < sim->nx+1; i += sim->nx) {
            printf("u[%d][%d] = %lf \t us = %lf --> Left and right\n", i, j, sim->U[i][j], sim->US[i][j]);
        }
    }

    // set identity for rectangle
    for (i = sim->n * D_IN; i < (D_IN + LBOX) * sim->n + 1; i++) {
        for (j = D_BOT * sim->n + 1; j < (D_BOT + 1) * sim->n + 1; j++) {
            printf("u[%d][%d] = %lf \t us = %lf --> Inside\n", i, j, sim->U[i][j], sim->US[i][j]);
        }
    }

    // set identity for upper and lower walls
    for (i = 0; i < sim->nx + 2; i++) {
        for (j = 0; j < sim->ny + 1; j += sim->ny) {
            printf("v[%d][%d] = %lf \t vs = %lf --> Upper an lower\n", i, j, sim->V[i][j], sim->VS[i][j]);
        }
    }

    // set identity for left and right walls
    for (j = 1; j < sim->ny; j++) {
        for (i = 0; i < sim->nx+2; i += sim->nx + 1) {
            printf("v[%d][%d] = %lf \t vs = %lf --> Left and right\n", i, j, sim->V[i][j], sim->VS[i][j]);
        }
    }

    // set identity for rectangle
    for (i = sim->n * D_IN + 1; i < (D_IN + LBOX) * sim->n + 1; i++) {
        for (j = D_BOT * sim->n; j < (D_BOT + 1) * sim->n + 1; j++) {
            printf("v[%d][%d] = %lf \t vs = %lf --> Inside\n", i, j, sim->V[i][j], sim->VS[i][j]);
        }
    }
}


void free_Sim_data(Sim_data *sim) {
    free(sim->u_data); free(sim->U);
    free(sim->v_data); free(sim->V);
    free(sim->p_data); free(sim->P);
#   if TEMP_MODE
    free(sim->T_data); free(sim->T);
#   endif
    free(sim);
}
