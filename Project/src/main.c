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
     RE = U0 Hbox / nu
    */

    printf("\n ============================ \e[1mLMECA2660 project in CFD\e[0m ============================\n");
#   if USE_ADI
    printf("\e[1m       Re = %.0lf    dx = %.3lf    dt = %.4lf    T simu = %.2lf   ADI enabled \e[0m\n\n", RE, sim->h, sim->dt, sim->tsim);
#   else
    printf("\e[1m       Re = %.0lf    dx = %.3lf    dt = %.4lf    T simu = %.2lf   ADI disabled \e[0m\n\n", RE, sim->h, sim->dt, sim->tsim);
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
    double exec_time = (double) (now - start_exec);
    double time_left = exec_time * (nt - t) / t;

    // int hours = exec_time / 3600;
    // int seconds = exec_time % 3600;
    // int minutes = seconds / 60;
    // seconds = seconds % 60;

    int seconds = ((int) exec_time) % 60;
    int minutes = ((int) exec_time) / 60;

    int s_left = ((int) time_left) % 60;
    int m_left = ((int) time_left) / 60;

    // printf("\r\e[93m\e[1m%3d%%\e[0m [%.*s%*s]   \e[91mPoisson %d its\e[0m    %*d/%*d   t_sim = %*.3lf s     \e[92m[%02dh%02d:%02d\e[0m",
    //        val, lpad, PBSTR, rpad, "", res, ndec1, t, ndec1, nt, ndec2+4, (double) t * dt, hours, minutes, seconds);

    printf("\r\e[93m\e[1m%3d%%\e[0m [%.*s%*s]   \e[91mPoisson %d its\e[0m    %*d/%*d   t_sim = %*.3lf s     \e[92m[%02d:%02d < %02d:%02d]\e[0m",
           val, lpad, PBSTR, rpad, "", res, ndec1, t, ndec1, nt, ndec2+4, (double) t * dt, minutes, seconds, m_left, s_left);

    fflush(stdout);
}


void integrate_flow(Sim_data *sim, Poisson_data *poisson, ADI_data *adi) {
    time_t start_exec;

    int nIterations;
    int t = 0;
    int ndec1 = (int) ceil(log10(sim->nt + 1));
    int ndec2 = (int) ceil(log10(sim->tsim + 1.));

    // set_bd_conditions(sim, sim->U, sim->V);
    // set_ghost_points(sim);
    
    save_fields(sim, t);
    time(&start_exec);

    for (t = 1; t <= sim->nt; t++) {
        
        // in what order should these three blocks be ???
        
        sim->t = t;
        set_mesh_velocity(sim);
        
        set_bd_conditions(sim, sim->U, sim->V);
        set_bd_conditions(sim, sim->US, sim->VS);
        set_ghost_points(sim);

        compute_convection(sim);

#       if USE_ADI
        predictor_step_u_adi(sim, adi);
        predictor_step_v_adi(sim, adi);
#       else
        predictor_step(sim);
#       endif

        nIterations = poisson_solver(sim, poisson);
        corrector_step(sim);

        // set_bd_conditions(sim, sim->U, sim->V);
        // set_ghost_points(sim);

        printProgress((double) t / (double) sim->nt, nIterations, t, sim->nt, sim->dt, ndec1, ndec2, start_exec);
        if ((t % SAVE_MODULO == 0) && (SAVE)) {
            save_fields(sim, t);
        }

        swap_next_previous(sim);
    }

    printf("\n\n");

}


int main(int argc, char *argv[]){
    // argv : ./cfd -ksp_type fgmres -pc_type lu -n 40 -dt 0.001 -tend 5.    

    int tmp1;
    int n = 40;
    double dt = 0.002;
    double tend = 10.;
    double tmp2;

    char *endptr;

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-n") == 0) {
            if (argc > i+1) tmp1 = (int) strtol(argv[i+1], &endptr, 10);
            if ((argc > i+1) && (*endptr == '\0')) {
                i++;
                n = tmp1;
            } else {
                printf("Could not read value for parameter [n]. It was set to 40 by default.\n");
                n = 40;
            }
        } else if (strcmp(argv[i], "-dt") == 0) {
            if (argc > i+1) tmp2 = (double) strtod(argv[i+1], &endptr);
            if ((argc > i+1) && (*endptr == '\0')) {
                i++;
                dt = tmp2;
            } else {
                printf("Could not read value for parameter [dt]. It was set to 0.002 by default.\n");
                dt = 0.002;
            }
        } else if (strcmp(argv[i], "-tend") == 0) {
            if (argc > i+1) tmp2 = (double) strtod(argv[i+1], &endptr);
            if ((argc > i+1) && (*endptr == '\0')) {
                i++;
                tend = tmp2;
            } else {
                printf("Could not read value for parameter [tend]. It was set to 10 by default.\n");
                tend = 10.;
            }
        } 
    }

    argc = 3;
    PetscInitialize(&argc, &argv, 0, 0);

    Sim_data *simulation = (Sim_data *)malloc(sizeof(Sim_data));
    init_Sim_data(simulation, n, dt, tend);
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
