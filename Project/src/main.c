#include "adi.h"
#include "poisson.h"
#include "project.h"

// Fix: implement adi for temperature // bof bof
// Improvement: implement streched grid, adaptive time
// Improvement: use input file instead of defines

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

    printf(
        "\n ==================================== \e[1mLMECA2660 project in CFD\e[0m "
        "====================================\n");
    const char *dt_status = (ADAPTIVE_DT == 1) ? "adaptive" : "fixed";
    const char *adi_status = (USE_ADI == 1) ? "enabled" : "disabled";
    const char *tmp_status = (TEMP_MODE > 0) ? "enabled" : "disabled";
    printf("\e[1m       Re = %.0lf    dx = %.3lf    T_simu = %.2lf   dt %s    Temp %s   ADI %s \e[0m \n\n", RE, sim->h,
           sim->tsim, dt_status, tmp_status, adi_status);
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

void printProgress(Sim_data *sim, int res, int ndec1, int ndec2, struct timespec start_exec) {
    double percentage = sim->tnow / sim->tsim;

    int val = (int)(percentage * 100);
    int lpad = (int)(percentage * PBWIDTH);
    int rpad = PBWIDTH - lpad;

    struct timespec clock_now;
    clock_gettime(CLOCK_REALTIME, &clock_now);

    double exec_time =
        (double)((1e9 * clock_now.tv_sec + clock_now.tv_nsec) - (1e9 * start_exec.tv_sec + start_exec.tv_nsec)) * 1.e-9;
    double time_left = exec_time * (sim->tsim / sim->tnow - 1);
    time_left = fmax(time_left, 0.);

    // int hours = exec_time / 3600;
    // int seconds = exec_time % 3600;
    // int minutes = seconds / 60;
    // seconds = seconds % 60;

    double seconds = fmod(exec_time, 60.);
    double minutes = floor(exec_time / 60.);
    double s_left = fmod(time_left, 60.);
    double m_left = floor(time_left / 60.);

    // printf("\r\e[93m\e[1m%3d%%\e[0m [%.*s%*s]   \e[91mPoisson %d its\e[0m    %*d/%*d   t_sim = %*.3lf s
    // \e[92m[%02dh%02d:%02d\e[0m",
    //        val, lpad, PBSTR, rpad, "", res, ndec1, t, ndec1, nt, ndec2+4, (double) t * dt, hours, minutes, seconds);

    char info_msg[250] =
        "\r\e[93m\e[1m%3d%%\e[0m [%.*s%*s]   \e[91mPoisson %d its\e[0m    t_sim = %*.3lf s    \e[34mdt = %.3lf ms < "
        "%.3lf ms\e[0m     Re_w,h = %4.1lf, %4.1lf    \e[92m[%02.0f:%02.0f < %02.0f:%02.0f] \e[0m ";
    printf(info_msg, val, lpad, PBSTR, rpad, "", res, ndec2 + 4, sim->tnow, 1.e3 * sim->dt, 1.e3 * sim->dt_stable,
           sim->rew, sim->reh, minutes, seconds, m_left, s_left);
    fflush(stdout);
}

void integrate_flow(Sim_data *sim, Poisson_data *poisson, ADI_data *adi) {
    struct timespec start_exec;

    int nIterations;
    int ndec1 = (int)ceil(log10(sim->nt + 1));
    int ndec2 = (int)ceil(log10(sim->tsim + 1.));
    sim->first_iteration = 1;

    // sim->t = 0;
    save_fields(sim);
    save_diagnostics(sim, 1);

    while (sim->tnow + 1e-9 < sim->tsim) {
        // CONVECTION
        compute_convection(sim);
        compute_convection_temperature(sim);  // does nothing if TEMP_MODE disabled

        // PREDICTOR STEP
#if USE_ADI
        predictor_step_u_adi(sim, adi);
        predictor_step_v_adi(sim, adi);
#else
        predictor_step(sim);
#endif

        // UPDATE TIME AND BOUNDARY CONDITIONS
        sim->tnow += sim->dt;
        sim->elapsed += sim->dt;

        set_mesh_velocity(sim, sim->tnow);
        set_bd_conditions(sim);
        set_ghost_points(sim);
        set_boundary_temperature(sim);  // does nothing if TEMP_MODE disabled

        // SOLVE POISSON EQUATION
        // clock_t start = clock();
        nIterations = poisson_solver(sim, poisson);
        // printf("Time to solve = %.3f s\n", (double) (clock() - start) / CLOCKS_PER_SEC);

        // CORRECTOR STEP
        corrector_step(sim);
        corrector_step_temperature(sim);  // does nothing if TEMP_MODE disabled

        // KEEP HISTORY OF HX, HY, HT
        swap_next_previous(sim);

        // PROGRESS BAR
        if (sim->first_iteration) {
            clock_gettime(CLOCK_REALTIME,
                          &start_exec);  // don't take first iteration into account, since much much longer
            sim->first_iteration = 0;
        }
        printProgress(sim, nIterations, ndec1, ndec2, start_exec);

        // SAVE RESULTS
        if ((sim->save_freq > 0.) && (sim->elapsed > sim->save_freq)) {
            save_fields(sim);
            save_diagnostics(sim, 1);
            sim->elapsed = fmod(sim->elapsed, sim->save_freq);
        } else {
            save_diagnostics(sim, 0);  // second argument is a flag to indicate if the fields are saved
        }
    }

    printf("\n\n");
}

int main(int argc, char *argv[]) {
    // ./cfd -ksp_type fgmres -pc_type lu -n 46 -dt 0.001 -tend 50. -freq 0.1 -dir new_case
    // -ksp_type : solver used by PETSC to solve the poisson equation of the two-step method
    // -pc_type  : pre-conditionner of the PETSC solver
    // -n        : number of spatial steps per distance H_box (must be >= 5)
    // -dt       : time step, automatically computed with CFL if set to 0.
    // -tend     : final time of the simulation
    // -freq     : save the fields u, v, p, T every ... "seconds", 0 to never save
    // -dir      : sub-directory in which the files are saved (inside the directory "myPath")

    /**
     * Handle the flags received in the command line
     */

    // Default values, but shouldn't be used in principle
    int n = 40;
    double dt = 0.002;
    double tend = 5.;
    double save_freq = 0.1;

    char myPath[50] = "./results/";
    char *endptr;

    if (argc != 15) {
        printf("Bad input  -->  usage: \n");
        printf("./cfd -ksp_type fgmres -pc_type lu -n 46 -dt 0.001 -tend 50. -freq 0.1 -dir new_case\n");
        return EXIT_FAILURE;
    }

    int i;
    int count_args = 0;
    for (i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-n") == 0) {
            if (argc > i + 1) n = (int)strtol(argv[i + 1], &endptr, 10);
            if ((argc > i + 1) && (*endptr == '\0') && (n >= 5)) {
                i++;
                count_args++;
            } else if (n < 5) {
                ABORT_MSG("Parameter [n] should be at least 5 !");
            } else {
                ABORT_MSG("Could not read value for parameter [n] !");
            }
        } else if (strcmp(argv[i], "-dt") == 0) {
            if (argc > i + 1) dt = (double)strtod(argv[i + 1], &endptr);
            if ((argc > i + 1) && (*endptr == '\0') && (dt > 0.)) {
                i++;
                count_args++;
            } else if (dt <= 0.) {
                ABORT_MSG("[dt] must be > 0 !");
                return EXIT_FAILURE;
            } else {
                ABORT_MSG("Could not read value for parameter [dt] !");
                return EXIT_FAILURE;
            }
        } else if (strcmp(argv[i], "-tend") == 0) {
            if (argc > i + 1) tend = (double)strtod(argv[i + 1], &endptr);
            if ((argc > i + 1) && (*endptr == '\0') && (tend > 0.)) {
                i++;
                count_args++;
            } else if (tend <= 0.) {
                ABORT_MSG("[tend] must be > 0. !");
                return EXIT_FAILURE;
            } else {
                ABORT_MSG("Could not read value for parameter [tend] !\n");
                return EXIT_FAILURE;
            }
        } else if (strcmp(argv[i], "-freq") == 0) {
            if (argc > i + 1) save_freq = (double)strtod(argv[i + 1], &endptr);
            if ((argc > i + 1) && (*endptr == '\0') && (save_freq >= 0.)) {
                i++;
                count_args++;
            } else if (save_freq < 0) {
                ABORT_MSG("[freq] must be > 0, or = 0 if you don't want to save the fields !");
                return EXIT_FAILURE;
            } else {
                ABORT_MSG("Could not read value for parameter [freq]");
                return EXIT_FAILURE;
            }
        } else if (strcmp(argv[i], "-dir") == 0) {
            if (argc > i + 1) {
                struct stat st = {0};
                char res;

                strcat(myPath, argv[i + 1]);
                strcat(myPath, "/");
                count_args++;

                if (stat(myPath, &st) == -1) {
                    mkdir(myPath, 0700);
                } else {
                    printf("The directory %s already exists. Do you want to overwrite it ? [y] or [n]: ", argv[i + 1]);
                    if ((scanf("%c", &res) == 1) && (res == 'y')) {
                        continue;
                    } else {
                        printf("Aborting the program\n");
                        return EXIT_FAILURE;
                    }
                }
            } else {
                printf("Could not read value for parameter [dir]. Writing files in parent directory.\n");
                return EXIT_FAILURE;
            }
        }
    }

    if (count_args != 5) {
        printf("Expected 5 options, only got %d  --->  Usage:\n", count_args);
        printf("./cfd -ksp_type fgmres -pc_type lu -n 46 -dt 0.001 -tend 50. -freq 0.1 -dir new_case\n");
        return EXIT_FAILURE;
    }

    /**
     * Initialization function of PETSC
     */
    int argc_petsc = 3;  // PETSC is not concerned by the other arguments
    PetscInitialize(&argc_petsc, &argv, 0, 0);

    /**
     * Initialize the simulation and the fields at t=0
     */
    Sim_data *simulation = (Sim_data *)malloc(sizeof(Sim_data));
    init_Sim_data(simulation, n, dt, tend, save_freq, myPath);
    init_fields(simulation);

    /**
     * Initialize the poisson solver and create the laplacian matrix of the system
     */
    // clock_t start = clock();
    Poisson_data *poisson = (Poisson_data *)malloc(sizeof(Poisson_data));
    initialize_poisson_solver(simulation, poisson);
    // printf("Time to create Matrix = %.3f s\n", (double) (clock() - start) / CLOCKS_PER_SEC);

    /**
     * Initialize the ADI solver if needed. Not sure boundary conditions are 100% correct yet.
     */
    ADI_data *adi_solver = (ADI_data *)malloc(sizeof(ADI_data));
    init_adi_solver(simulation, adi_solver);  // does nothing if ADI is disabled

    /**
     * MAIN PROCESS
     */
#if TEST_POISSON
    test_poisson(simulation, poisson);
#else
    display_info(simulation, "full");  // "full"
    integrate_flow(simulation, poisson, adi_solver);
#endif

    /**
     * FREE THE MEMORY
     */
    free_Sim_data(simulation);
    free_poisson_solver(poisson);
    free_adi_solver(adi_solver);
    PetscFinalize();

    return EXIT_SUCCESS;
}
