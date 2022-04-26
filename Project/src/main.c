#include "project.h"
#include "poisson.h"
#include "adi.h"

// TODO: implement adi for temperature // bof bof
// TODO: implement streched grid

void display_info(Sim_data *sim, char *mode) {
    /*
    Choice of Nondimentionalization:
     x_ = x / Hbox
     y_ = y / Hbox
     u_ = u / U0
     v_ = v / U0
     t_ = t / (Hbox/U0)
     p_ = (p-p_ref) / (rho U0^2)
     T_ = (T - T0) / (T1 - T0)
     RE = U0 Hbox / nu
    */

    printf("\n =================================== \e[1mLMECA2660 project in CFD\e[0m ===================================\n");
    const char *adi_status = (USE_ADI == 1) ? "enabled" : "disabled";
    const char *tmp_status = (TEMP_MODE > 0) ? "enabled" : "disabled";
    printf("\e[1m       Re = %.0lf    dx = %.3lf    dt = %.4lf    T_simu = %.2lf   Temp %s   ADI %s \e[0m \n\n",
           RE, sim->h, sim->dt, sim->tsim, tmp_status, adi_status);
    char description[] = {
        "        ╭──────────────────────────────────────────────────────────────╮\n"
        "        │                       %2d                                     │\n"
        "        │              <------------------->                           │\n"
        "        │                                                              │\n"
        "        │     %2d      ╭─────────────────────╮            %2d            │\n"
        "        │  <------->  │              |%2d    │  <---------------------> │\n"
        "        │             ╰─────────────────────╯                          │\n"
        "        │                        |                                     │\n"
        "        │                        |%2d                                   │\n"
        "        │                        |                                     │\n"
        "        ╰──────────────────────────────────────────────────────────────╯\n\n"};
    
    if (strcmp(mode, "full") == 0) {
        printf(description, LBOX, D_IN, L_ - D_IN - LBOX, 1, D_BOT);
    }
}


void printProgress(double percentage, int res, int t, int nt, double dt, int ndec1, int ndec2, struct timespec start_exec) {
    int val = (int) (percentage * 100);
    int lpad = (int) (percentage * PBWIDTH);
    int rpad = PBWIDTH - lpad;

    struct timespec now;
    clock_gettime(CLOCK_REALTIME, &now);

    double exec_time = (double) ((1e9 * now.tv_sec + now.tv_nsec) - (1e9 * start_exec.tv_sec + start_exec.tv_nsec)) * 1.e-9;
    double time_left = exec_time * (nt - t) / t;

    // int hours = exec_time / 3600;
    // int seconds = exec_time % 3600;
    // int minutes = seconds / 60;
    // seconds = seconds % 60;

    double seconds = fmod(exec_time, 60.);
    double minutes = exec_time / 60.;
    double s_left = fmod(time_left, 60.);
    double m_left = time_left / 60.;

    // printf("\r\e[93m\e[1m%3d%%\e[0m [%.*s%*s]   \e[91mPoisson %d its\e[0m    %*d/%*d   t_sim = %*.3lf s     \e[92m[%02dh%02d:%02d\e[0m",
    //        val, lpad, PBSTR, rpad, "", res, ndec1, t, ndec1, nt, ndec2+4, (double) t * dt, hours, minutes, seconds);

    printf("\r\e[93m\e[1m%3d%%\e[0m [%.*s%*s]   \e[91mPoisson %d its\e[0m    %*d/%*d   t_sim = %*.3lf s     \e[92m[%02.0f:%02.0f < %02.0f:%02.0f] \e[0m ",
           val, lpad, PBSTR, rpad, "", res, ndec1, t, ndec1, nt, ndec2+4, (double) t * dt, minutes, seconds, m_left, s_left);

    fflush(stdout);
}


void integrate_flow(Sim_data *sim, Poisson_data *poisson, ADI_data *adi) {
    struct timespec start_exec;
    
    int nIterations;
    int ndec1 = (int) ceil(log10(sim->nt + 1));
    int ndec2 = (int) ceil(log10(sim->tsim + 1.));
    
    sim->t = 0;    
    save_fields(sim, sim->t);
    write_diagnostics(sim, sim->t);

    while (sim->t < sim->nt) {

        // CONVECTION        
        compute_convection(sim);
#       if TEMP_MODE
        compute_convection_temperature(sim);
#       endif
        
        // PREDICTOR STEP
#       if USE_ADI
        predictor_step_u_adi(sim, adi);
        predictor_step_v_adi(sim, adi);
#       else
        predictor_step(sim);
#       endif

        // UPDATE TIME AND BOUNDARY CONDITIONS
        sim->t = sim->t + 1;
        set_mesh_velocity(sim, sim->t * sim->dt);
        set_bd_conditions(sim, sim->U, sim->V);   // sets u and v at t=(n+1)  # no influence on corrector step
        set_bd_conditions(sim, sim->US, sim->VS);
        set_ghost_points(sim);
#       if TEMP_MODE
        set_boundary_temperature(sim);
#       endif

        // SOLVE POISSON EQUATION
        // clock_t start = clock();
        nIterations = poisson_solver(sim, poisson);
        // printf("Time to solve = %.3f s\n", (double) (clock() - start) / CLOCKS_PER_SEC);

        // CORRECTOR STEP
        corrector_step(sim);
#       if TEMP_MODE
        corrector_step_temperature(sim);
#       endif

        // KEEP HISTORY OF HX, HY, HT
        swap_next_previous(sim);

        // SAVE RESULTS
        if (sim->t == 1) clock_gettime(CLOCK_REALTIME, &start_exec);  // don't take first iteration into account, since much much longer
        printProgress((double) sim->t / (double) sim->nt, nIterations, sim->t, sim->nt, sim->dt, ndec1, ndec2, start_exec);
        write_diagnostics(sim, sim->t);
        if ((sim->t % SAVE_MODULO == 0) && (SAVE)) {
            save_fields(sim, sim->t);
        }
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

    int argc_petsc = 3;
    PetscInitialize(&argc_petsc, &argv, 0, 0);

    Sim_data *simulation = (Sim_data *)malloc(sizeof(Sim_data));
    init_Sim_data(simulation, n, dt, tend);
    init_fields(simulation);

    clock_t start = clock();
    Poisson_data *poisson = (Poisson_data *)malloc(sizeof(Poisson_data));
    initialize_poisson_solver(simulation, poisson);
    printf("Time to create Matrix = %.3f s\n", (double) (clock() - start) / CLOCKS_PER_SEC);

    ADI_data *adi_solver = (ADI_data *)malloc(sizeof(ADI_data));
#   if USE_ADI
    init_adi_solver(simulation, adi_solver);
#   endif


    // MAIN PROCESS
#   if TEST_POISSON
    test_poisson(simulation, poisson);
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
