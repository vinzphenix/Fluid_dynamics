#include "thomas.h"
#include "fd_solver.c"


int init_data_sim(data_Sim *sim) {

#   if SAVE == 1
        sprintf(filename, "%ssolution.txt", path);
#   elif (SAVE == 2) && (WAVEPACKET == 0)
        sprintf(filename, "%ssolution_%c%c_%d.txt", path, SCHEME_A, SCHEME_B, N);
#   elif (SAVE == 2) && (WAVEPACKET == 1)
        sprintf(filename, "%swavepacket_%c%c_%d.txt", path, SCHEME_A, SCHEME_B, N);
#   elif (SAVE == 3) && (WAVEPACKET == 0)
        sprintf(filename, "%snonuniform_%c%c_%d.txt", path, SCHEME_A, SCHEME_B, N);
#   endif

    sim->h = L / N;
    sim->dt = CFL * sim->h * (1 - A) / fabs(C);
    sim->M = ceil(TEND / sim->dt);

#   if (SCHEME_A == 'I')
        sim->u = (double *)calloc(8 * N, sizeof(double));
        sim->x1 = sim->u + 5 * N;
        sim->at = sim->u + 6 * N;
        sim->q = sim->u + 7 * N;
#   else
        sim->u = (double *)calloc(5 * N, sizeof(double));
#   endif
    // malloc would work as good as calloc since we never use these arrays uninitialized 
    // except "du", but it is multiplied by 0 in that case

    sim->ul = sim->u + N;
    sim->us = sim->u + 2 * N;
    sim->du = sim->u + 3 * N;

    sim->dg = sim->u + 4 * N;
    double xi;
    for (int i = 0; i < N; i++) {
        xi = -L / 2. + i * sim->h;
        sim->dg[i] = 1. - A * cos(2 * M_PI * xi / L);
    }

    return 0;
}

void free_data_sim(data_Sim *sim) {
    free(sim->u);
    free(sim);
}

void save_array(data_Sim *sim, int t) {
    FILE *ptr;
    if (t == 0) {
        ptr = fopen(filename, "w");
        fprintf(ptr, "%c%c\n", SCHEME_A, SCHEME_B);
        fprintf(ptr, "%lf %lf %lf %lf %lf %lf\n", C, SIGMA, UMAX, L, sim->dt, A);
    } else {
        ptr = fopen(filename, "a");
    }

    fprintf(ptr, "%.5lf ", t * sim->dt);
    for (int i = 0; i < N; i++) {
        fprintf(ptr, "%.5e ", sim->u[i]);
    }

    fprintf(ptr, "\n");
    fclose(ptr);
}

void set_u_initial(data_Sim *sim) {
    double x;
    for (int i = 0; i < N; i++) {
        x = -L / 2. + i * sim->h;
        x = x - A * L / (2 * M_PI) * sin(2 * M_PI * x / L);
        sim->u[i] = UMAX * exp(-x * x / (SIGMA * SIGMA)) * sim->dg[i];
#       if WAVEPACKET > 0
        sim->u[i] *= cos(2 * M_PI * 16 * x / L);
#       endif
    }

#   if SAVE
        save_array(sim, 0);
#   endif
}

void display_diagnostic(data_Sim *sim, int t_idx) {
    double I = 0.;
    double E = 0.;
    double R = 0.;
    double x, u_exact, arg;
    double t = sim->dt * t_idx;

    for (int i = 0; i < N; i++) {
        x = -L / 2. + i * sim->h;  // uniform
        x = x - A * L / (2 * M_PI) * sin(2 * M_PI * x / L);  // scaling
        arg = fmod(x - C * t - L / 2., L) + L / 2.;  // working for C > 0
        u_exact = UMAX * exp(-pow(arg / SIGMA, 2.)) * sim->dg[i];

#       if WAVEPACKET > 0
        u_exact *= cos(2 * M_PI * 16 * (x - C * t) / L);
#       endif

        I += sim->u[i];
        E += pow(sim->u[i], 2.);
        R += pow(sim->u[i] - u_exact, 2.);
    }

    I *= sim->h / (SIGMA * UMAX);
    E *= sim->h / (SIGMA * UMAX * UMAX);
    R *= sim->h / (SIGMA * UMAX * UMAX);
    printf("iteration %3d   t = %.3f  : \t  I = %-7.3lf  E = %-7.3lf  R = %-7.3f \n", t_idx, t, I, E, R);
}

void RK4C(data_Sim *sim) {
    double *u = sim->u;
    double *us = sim->us;
    double *ul = sim->ul;
    double *du = sim->du;

    int i, k, t;

    for (t = 1; t <= sim->M; t++) {

        memcpy(us, u, N * sizeof(double));
        for (k = 0; k < 4; k++) {
            for (i = 0; i < N; i++) {
                ul[i] = (us[i] + ALPHA[k] * sim->dt * du[i]) / sim->dg[i];
            }
            f_eval(sim);
            for (i = 0; i < N; i++) {
                u[i] += GAMMA[k] * sim->dt * du[i];
            }
        }

#       if SAVE > 0
            display_diagnostic(sim, t);
            save_array(sim, t);
#       endif
        
    }

    return;
}

int main(int argc, char *argv[]) {

    data_Sim *simulation = (data_Sim *)malloc(sizeof(data_Sim));
    init_data_sim(simulation);
    set_u_initial(simulation);

    clock_t start = clock();
    RK4C(simulation);
    clock_t end = clock();
    printf("time taken = %.3lf ms\n", 1000 * ((double)end - (double)start) / CLOCKS_PER_SEC);

    free_data_sim(simulation);
    return 0;
}

// gcc -o hw main.c -lm -O3 && ./hw
