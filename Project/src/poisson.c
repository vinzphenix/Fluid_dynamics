
#include <mpi.h>
#include "poisson.h"

/*Called by poisson_solver at each time step*/
void computeRHS(Sim_data *sim, double *rhs, PetscInt rowStart, PetscInt rowEnd) {

    int i, j;
    int r = rowStart;

    /*int i_w_left =  D_IN * sim->n - 1;
    int i_w_right = (D_IN + LBOX) * sim->n;
    int j_w_below = D_BOT * sim->n - 1;
    int j_w_above = (D_BOT + 1) * sim->n;*/
    
#   if NO_SLIP
    // YOU MUST IMPOSE A ZERO-MASS FLOW HERE ...
    i = sim->nx;
    double flow_in = 0.;
    double flow_out = 0.;
    for (j = 1; j < sim->ny + 1; j++) {
        flow_in += sim->US[0][j];
        flow_out += sim->US[i][j];
    }
    // if (sim->elapsed > sim->save_freq)
    //     printf("\nBEFORE : In : %15.8lf  Out : %15.8lf\n", flow_in, flow_out);

    for (j = 1; j < sim->ny + 1; j++) {
        sim->US[i][j] = sim->US[i][j] - sim->uin[j] * (flow_out - flow_in) / sim->ny; 
        sim->U[i][j] = sim->US[i][j];
    }
#   endif    
    
    // Going through the whole domain, included the points inside, because not taken into account later anyway
    for (i = 0; i < sim->nx; i++) {
        for (j = 0; j < sim->ny; j++) {
            // 1/h factor of the divergence operator is taken into account in the matrix
            rhs[r++] = ((sim->US[i+1][j+1] - sim->US[i  ][j+1]) + (sim->VS[i+1][j+1] - sim->VS[i+1][j  ])) / sim->dt;
        }
    }

    /*Do not forget that the solution for the Poisson equation is defined within a constant.
    One point from Phi must then be set to an abritrary constant.*/
    rhs[sim->idx_set_p] = 0.;  // set value of phi at inflow
}


void computeRHS_TEST(Sim_data *sim, double *rhs, PetscInt rowStart, PetscInt rowEnd) {
    int i, j;
    int r = rowStart;

    // be carefull : (integral_boundary fluxes) MUST BE EQUAL TO (int_domain) source term (here 0.)
    for (i = 0; i < sim->nx; i++) {
        for (j = 0; j < sim->ny; j++, r++) {
            if ((j == 0)) {
                rhs[r] += 2. / sim->h * cos(2. * M_PI * (i + 0.5) * sim->h / L_);
            }
            if ((j == sim->ny-1)) {
                rhs[r] += -2. / sim->h * sin(4. * M_PI * (i + 0.5) * sim->h / L_);
            }
        }
    }
    rhs[0] = 0.;
}


/*To call at each time step after computation of U_star. This function solves the poisson equation*/
/*and copies the solution of the equation into your vector Phi*/
int poisson_solver(Sim_data *sim, Poisson_data *data) {

    /* Solve the linear system Ax = b for a 2-D poisson equation on a structured grid */
    int its;
    PetscInt rowStart, rowEnd;
    PetscScalar *rhs, *sol;

    KSP sles = data->sles;
    Vec b = data->b;
    Vec x = data->x;

    /* Fill the right-hand-side vector : b */
    VecGetOwnershipRange(b, &rowStart, &rowEnd);
    VecGetArray(b, &rhs);
#   if TEST_POISSON
    computeRHS_TEST(sim, rhs, rowStart, rowEnd); /*MODIFY THE PROTOTYPE HERE*/
#   else
    computeRHS(sim, rhs, rowStart, rowEnd); /*MODIFY THE PROTOTYPE HERE*/
#   endif
    VecRestoreArray(b, &rhs);


    /*Solve the linear system of equations */
    KSPSolve(sles, b, x);  // printf("res = %d\n", res);
    KSPGetIterationNumber(sles, &its);
    // if (sim->t % SAVE_MODULO == 0)
    // printf("t = %5d \t Poisson solved in %d iterations \r", sim->t, its);
    VecGetArray(x, &sol);

    int q = sim->size_p;
    int r;
    for(r = rowStart; r < rowEnd; r++, q++){
        sim->p_data[q] = sol[r];
    }

    VecRestoreArray(x, &sol);
    free(sol);
    free(rhs);
    
    return its;
}

/*
This function is called only once during the simulation, i.e. in initialize_poisson_solver.
*/
void computeLaplacianMatrixNOIF(Sim_data *sim, Mat A, int rowStart, int rowEnd) {
    int i, j, idx, block;
    int k = sim->ny;
    double alpha = sim->n;  // 1/h^2  *  h  =  n
    
    int i_w_left =  D_IN * sim->n - 1;
    int i_w_right = (D_IN + LBOX) * sim->n;
    int j_w_below = D_BOT * sim->n - 1;
    int j_w_above = (D_BOT + 1) * sim->n;

    int i_start[4] = {1, i_w_left, i_w_left, i_w_right + 1}; // included
    int i_final[4] = {i_w_left, i_w_right + 1, i_w_right + 1, sim->nx - 1}; // not included
    int j_start[4] = {1, 1, j_w_above + 1, 1};
    int j_final[4] = {sim->ny - 1, j_w_below, sim->ny - 1, sim->ny - 1};

    // main blocks no touching boundaries
    for (block = 0; block < 4; block++) {
        for (i = i_start[block]; i < i_final[block]; i++) {
            for (j = j_start[block]; j < j_final[block]; j++) {
                idx = i * k + j;
                MatSetValue(A, idx, idx+1, alpha, INSERT_VALUES);
                MatSetValue(A, idx, idx-1, alpha, INSERT_VALUES);
                MatSetValue(A, idx, idx+k, alpha, INSERT_VALUES);
                MatSetValue(A, idx, idx-k, alpha, INSERT_VALUES);
                MatSetValue(A, idx, idx, -4. * alpha, INSERT_VALUES);
            }
        }
    }

    // external corners of rectangle have also 4 neighbours (using the classic stencil)
    for (i = i_w_left; i < i_w_right + 1; i += (i_w_right - i_w_left)) {
        for (j = j_w_below; j < j_w_above + 1; j += (j_w_above - j_w_below)) {
            idx = i * k + j;
            MatSetValue(A, idx, idx+1, alpha, INSERT_VALUES);
            MatSetValue(A, idx, idx-1, alpha, INSERT_VALUES);
            MatSetValue(A, idx, idx+k, alpha, INSERT_VALUES);
            MatSetValue(A, idx, idx-k, alpha, INSERT_VALUES);
            MatSetValue(A, idx, idx, -4. * alpha, INSERT_VALUES);
        }
    }

    // inside rectangle
    for (i = i_w_left + 1; i < i_w_right; i++) {
        for (j = j_w_below + 1; j < j_w_above; j++) {
            idx = i * k + j;
            MatSetValue(A, idx, idx, alpha, INSERT_VALUES);
        }
    }

    // only one neighbour missing above or below
    i_start[0] = 1; i_final[0] = sim->nx - 1;          // loop index : start
    i_start[1] = i_w_left + 1; i_final[1] = i_w_right; // loop index : end
    j_start[0] = 0; j_final[0] = sim->ny - 1;          // one shot
    j_start[1] = j_w_above; j_final[1] = j_w_below;    // one shot

    for (block = 0; block < 2; block++) {      
        for (i = i_start[block]; i < i_final[block]; i++) {

            j = j_start[block]; // below missing
            idx = i * k + j;
            MatSetValue(A, idx, idx+1, alpha, INSERT_VALUES);
            MatSetValue(A, idx, idx+k, alpha, INSERT_VALUES);
            MatSetValue(A, idx, idx-k, alpha, INSERT_VALUES);
            MatSetValue(A, idx, idx, -3. * alpha, INSERT_VALUES);

            j = j_final[block]; // above missing
            idx = i * k + j;
            MatSetValue(A, idx, idx-1, alpha, INSERT_VALUES);
            MatSetValue(A, idx, idx+k, alpha, INSERT_VALUES);
            MatSetValue(A, idx, idx-k, alpha, INSERT_VALUES);
            MatSetValue(A, idx, idx, -3. * alpha, INSERT_VALUES);

        }
    }

    // only one neighbour missing left or right
    j_start[0] = 1; j_final[0] = sim->ny - 1;           // loop index
    j_start[1] = j_w_below + 1; j_final[1] = j_w_above; // loop index
    i_start[0] = 0; i_final[0] = sim->nx - 1;           // one shot
    i_start[1] = i_w_right; i_final[1] = i_w_left;     // one shot

    for (block = 0; block < 2; block++) {      
        
        i = i_start[block]; // left missing
        for (j = j_start[block]; j < j_final[block]; j++) {
            idx = i * k + j;
            MatSetValue(A, idx, idx-1, alpha, INSERT_VALUES);
            MatSetValue(A, idx, idx+1, alpha, INSERT_VALUES);
            MatSetValue(A, idx, idx+k, alpha, INSERT_VALUES);
            MatSetValue(A, idx, idx, -3. * alpha, INSERT_VALUES);
        }

        i = i_final[block]; // right missing
        for (j = j_start[block]; j < j_final[block]; j++) {
            idx = i * k + j;
            MatSetValue(A, idx, idx-1, alpha, INSERT_VALUES);
            MatSetValue(A, idx, idx+1, alpha, INSERT_VALUES);
            MatSetValue(A, idx, idx-k, alpha, INSERT_VALUES);
            MatSetValue(A, idx, idx, -3. * alpha, INSERT_VALUES);
        }
    }

    // lower left corner (used to set reference)
    // idx = 0;
    // MatSetValue(A, idx, idx, alpha, INSERT_VALUES);

    // set reference
    idx = sim->idx_set_p;
    MatSetValue(A, idx, idx-1, 0., INSERT_VALUES);
    MatSetValue(A, idx, idx+1, 0., INSERT_VALUES);
    MatSetValue(A, idx, idx+k, 0., INSERT_VALUES);
    MatSetValue(A, idx, idx, alpha, INSERT_VALUES);

    // lower left
    idx = 0 * k + 0;
    MatSetValue(A, idx, idx, -2. * alpha, INSERT_VALUES);
    MatSetValue(A, idx, idx+1, alpha, INSERT_VALUES);
    MatSetValue(A, idx, idx+k, alpha, INSERT_VALUES);

    // upper left corner
    idx = 0 * k + (k-1);
    MatSetValue(A, idx, idx, -2. * alpha, INSERT_VALUES);
    MatSetValue(A, idx, idx-1, alpha, INSERT_VALUES);
    MatSetValue(A, idx, idx+k, alpha, INSERT_VALUES);

    // lower right corner
    idx = (sim->nx - 1) * k + 0;
    MatSetValue(A, idx, idx, -2. * alpha, INSERT_VALUES);
    MatSetValue(A, idx, idx+1, alpha, INSERT_VALUES);
    MatSetValue(A, idx, idx-k, alpha, INSERT_VALUES);

    // upper right corner
    idx = (sim->nx - 1) * k + (k-1);
    MatSetValue(A, idx, idx, -2. * alpha, INSERT_VALUES);
    MatSetValue(A, idx, idx-1, alpha, INSERT_VALUES);
    MatSetValue(A, idx, idx-k, alpha, INSERT_VALUES);
}

void computeLaplacianMatrix(Sim_data *sim, Mat A, int rowStart, int rowEnd) {

    int i, j, idx;
    int flag_right, flag_left, flag_above, flag_below;
    int k = sim->ny;
    double diag_value;
    double alpha = sim->n;  // (1/h^2)  *  h  =  n
    
    int i_w_left =  D_IN * sim->n - 1;
    int i_w_right = (D_IN + LBOX) * sim->n;
    int j_w_below = D_BOT * sim->n - 1;
    int j_w_above = (D_BOT + 1) * sim->n;

    for (i = 0; i < sim->nx; i++) {
        
        for (j = 0; j < sim->ny; j++) {
            idx = i * k + j;
            diag_value = 0.;

            flag_left  = (i == i_w_left)  && (j_w_below < j) && (j < j_w_above);
            flag_right = (i == i_w_right) && (j_w_below < j) && (j < j_w_above);
            flag_below = (j == j_w_below) && (i_w_left  < i) && (i < i_w_right);
            flag_above = (j == j_w_above) && (i_w_left  < i) && (i < i_w_right);

            if ((i_w_left < i) && (i < i_w_right) && (j_w_below < j) && (j < j_w_above)) { // inside rectangle
                diag_value += alpha;  // could use 1., but matrix conditionning would be worse, I supose
            } else { // outside rectangle
                if ((i != 0) && (!flag_right)) { // it has left neighbor
                    MatSetValue(A, idx, idx-k, alpha, INSERT_VALUES);
                    diag_value -= alpha;
                }
                if ((i != sim->nx - 1) && (!flag_left)) { // it has right neighbor
                    MatSetValue(A, idx, idx+k, alpha, INSERT_VALUES);
                    diag_value -= alpha;
                }
                if ((j != 0) && (!flag_above)) { // it has neighbor below
                    MatSetValue(A, idx, idx-1, alpha, INSERT_VALUES);
                    diag_value -= alpha;                }
                if ((j != sim->ny - 1) && (!flag_below)) { // it has neighbor above
                    MatSetValue(A, idx, idx+1, alpha, INSERT_VALUES);
                    diag_value -= alpha;
                }
            }
            MatSetValue(A, idx, idx, diag_value, INSERT_VALUES);
        }
    }

    /*Be careful; the solution from the system solved is defined within a constant.
    One point from Phi must then be set to an abritrary constant.*/

    idx = sim->idx_set_p;
    
    MatSetValue(A, idx, idx+k, 0., INSERT_VALUES);
    MatSetValue(A, idx, idx, alpha, INSERT_VALUES);

    if (0 < idx)  // not below
        MatSetValue(A, idx, idx - 1, 0., INSERT_VALUES);
    if (idx < sim->ny-1) // not above
        MatSetValue(A, idx, idx + 1, 0., INSERT_VALUES);
}


/*To call during the initialization of your solver, before the begin of the time loop*/
PetscErrorCode initialize_poisson_solver(Sim_data *sim, Poisson_data *data) {
    PetscInt rowStart = 0;         /*rowStart = 0*/
    PetscInt rowEnd = sim->size_p; /*rowEnd = the number of unknows*/
    PetscErrorCode ierr;

    int nphi = sim->size_p; /*WRITE HERE THE NUMBER OF UNKNOWS*/

    /* Create the right-hand-side vector : b */
    VecCreate(PETSC_COMM_WORLD, &(data->b));
    VecSetSizes(data->b, PETSC_DECIDE, nphi);
    VecSetType(data->b, VECSTANDARD);

    /* Create the solution vector : x */
    VecCreate(PETSC_COMM_WORLD, &(data->x));
    VecSetSizes(data->x, PETSC_DECIDE, nphi);
    VecSetType(data->x, VECSTANDARD);

    /* Create and assemble the Laplacian matrix : A  */
    MatCreate(PETSC_COMM_WORLD, &(data->A));
    MatSetSizes(data->A, PETSC_DECIDE, PETSC_DECIDE, nphi, nphi);
    MatSetType(data->A, MATAIJ);
    MatSeqAIJSetPreallocation(data->A, 5, NULL); // 5 nnz per row seems ok ! /*SET HERE THE NUMBER OF NON-ZERO DIAGONALS*/
    MatGetOwnershipRange(data->A, &rowStart, &rowEnd);

    // clock_t start = clock();
    // computeLaplacianMatrix_NO_IF(sim, data->A, rowStart, rowEnd); // not even faster
    computeLaplacianMatrix(sim, data->A, rowStart, rowEnd);
    // clock_t end = clock();

    ierr = MatAssemblyBegin(data->A, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);
    ierr = MatAssemblyEnd(data->A, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);

    /* Create the Krylov context */
    KSPCreate(PETSC_COMM_WORLD, &(data->sles));
    KSPSetOperators(data->sles, data->A, data->A);
    KSPSetType(data->sles, KSPFGMRES); // KSPGMRES seems the best, it will not be used if PC LU.
    PC prec;
    KSPGetPC(data->sles, &prec);
    PCSetType(prec, PCLU);
    KSPSetFromOptions(data->sles); // to uncomment if we want to specify the solver to use in command line. Ex: mpirun -ksp_type gmres -pc_type gamg
    KSPSetTolerances(data->sles, 1.e-12, 1e-12, PETSC_DEFAULT, PETSC_DEFAULT);
    KSPSetReusePreconditioner(data->sles, PETSC_TRUE);
    KSPSetUseFischerGuess(data->sles, 1, 4);
    KSPGMRESSetPreAllocateVectors(data->sles);

    // printf("Assembly of Mattrix and Vectors is done (time taken = %.3lf ms) \n\n", 1000 * ((double)end - (double)start) / CLOCKS_PER_SEC);

    return ierr;
}

/*To call after the simulation to free the vectors needed for Poisson equation*/
void free_poisson_solver(Poisson_data* data) {
    MatDestroy(&(data->A));
    VecDestroy(&(data->b));
    VecDestroy(&(data->x));
    KSPDestroy(&(data->sles));
    free(data);
}


void test_poisson(Sim_data *sim, Poisson_data *poisson) {
    int i, j;
    poisson_solver(sim, poisson);
    FILE *ptr = fopen("./results/test_poisson.txt", "w");
    fprintf(ptr, "%d\n", sim->n);
    for (i = 0; i < sim->nx; i++) {
        for (j = 0; j < sim->ny; j++) {
            fprintf(ptr, "%.4le ", sim->PHI[i][j]);
        }
        fprintf(ptr, "\n");
    } 
    fclose(ptr);
}
