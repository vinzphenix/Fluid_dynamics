#include "project.h"
#include "poisson.h"
#include "adi.h"


void display_info(Sim_data *sim, char *mode) {
    /*
    Choice of Nondimentionalization:
     x_ = x / Hbox
     y_ = y / Hbox
     u_ = u / U0
     v_ = v / U0
     t_ = t / (Hbox/U0)
     p_ = (p-p_ref) / (rho U0^2)
    */

    printf("\n ============================ \e[1mLMECA2660 project in CFD\e[0m ============================\n");
#   if USE_ADI
    printf("\e[1m       Re = %.0lf    dx = %.3lf    dt = %.4lf    T simu = %.2lf   ADI enabled \e[0m\n\n", RE, sim->h, sim->dt, TSIM);
#   else
    printf("\e[1m       Re = %.0lf    dx = %.3lf    dt = %.4lf    T simu = %.2lf   ADI disabled \e[0m\n\n", RE, sim->h, sim->dt, TSIM);
#   endif
    char description[] = {
        "        ╭──────────────────────────────────────────────────────────────╮\n"
        "        │                        %d                                     │\n"
        "        │              <------------------->                           │\n"
        "        │                                                              │\n"
        "        │      %d      ╭─────────────────────╮             %d            │\n"
        "        │  <------->  │              | %d    │  <---------------------> │\n"
        "        │             ╰─────────────────────╯                          │\n"
        "        │                        |                                     │\n"
        "        │                        | %d                                   │\n"
        "        │                        |                                     │\n"
        "        ╰──────────────────────────────────────────────────────────────╯\n\n"};
    
    if (strcmp(mode, "full") == 0) {
        printf(description, LBOX, D_IN, L_ - D_IN - LBOX, 1, D_BOT);
    }
}


void printProgress(double percentage, int res, int t, int nt, double dt, int ndec1, int ndec2, time_t start_exec) {
    int val = (int) (percentage * 100);
    int lpad = (int) (percentage * PBWIDTH);
    int rpad = PBWIDTH - lpad;

    time_t now;
    time(&now);
    int exec_time = (int) (now - start_exec);
    int hours = exec_time / 3600;
    int seconds = exec_time % 3600;
    int minutes = seconds / 60;
    seconds = seconds % 60;

    // int hours = 1;
    // int minutes = 1;
    // int seconds = 1;

    printf("\r\e[93m\e[1m%3d%%\e[0m [%.*s%*s]   \e[91mPoisson %d its\e[0m    %*d/%*d   t_sim = %*.3lf s     \e[92mexec: %02dh%02d:%02d\e[0m",
           val, lpad, PBSTR, rpad, "", res, ndec1, t, ndec1, nt, ndec2+4, (double) t * dt, hours, minutes, seconds);
    fflush(stdout);
}


void integrate_flow(Sim_data *sim, Poisson_data *poisson, ADI_data *adi) {
    time_t start_exec;

    int res;
    int t = 0;
    int ndec1 = (int) ceil(log10(sim->nt + 1));
    int ndec2 = (int) ceil(log10(TSIM + 1.));

    set_bd_conditions(sim, sim->U, sim->V);
    set_ghost_points(sim);
    save_fields(sim, t);

    time(&start_exec);

    for (t = 1; t <= sim->nt; t++) {
        
        sim->t = t;
        if (SIWNG_START < t * sim->dt) {
            sim->uMesh = ALPHA * sin(2. * M_PI * STROUHAL * (t * sim->dt - SIWNG_START));
        } else {
            sim->uMesh = 0.;
        }
        if ((PERT_START <= t * sim->dt) && (t * sim->dt <= PERT_START + PERT_DT)) {
            sim->vMesh = 0.5 * sin(2. * M_PI * (t * sim->dt - PERT_START) / PERT_DT);
        } else {
            sim->vMesh = 0.;
        }
        
        compute_convection(sim);

#       if USE_ADI
        predictor_step_u_adi(sim, adi);
        predictor_step_v_adi(sim, adi);
#       else
        predictor_step(sim);
#       endif

        res = poisson_solver(sim, poisson);
        corrector_step(sim);

        set_bd_conditions(sim, sim->U, sim->V);
        set_ghost_points(sim);

        swap_next_previous(sim);

        printProgress((double) t / (double) sim->nt, res, t, sim->nt, sim->dt, ndec1, ndec2, start_exec);
        if ((t % SAVE_MODULO == 0) && (SAVE)) {
            save_fields(sim, t);
        }
    }

    printf("\n\n");

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
#   elif TEST_TRIDIAGONAL
    double a[5] = {0.0, 1.0,-1.1, 1.2,-1.3};
    double b[5] = {5.0, 6.0, 7.0, 8.0, 9.0};
    double c[5] = {1.0,-1.1, 1.2,-1.3, 0.0};
    double q[5] = {1.0, 1.0, 1.0, 1.0, 1.0};
    solve_thomas(5, a, b, c, q);
    for (int i = 0; i < 5; i++) printf("q[%d] = %lf\n", i, q[i]);
#   else
    display_info(simulation, "no");
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