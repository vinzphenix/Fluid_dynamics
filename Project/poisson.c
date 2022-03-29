
#include <mpi.h>
#include "poisson.h"

/*Called by poisson_solver at each time step*/
/*More than probably, you should need to add arguments to the prototype ... */
/*Modification to do :*/
/*    -Impose zero mass flow here by changing value of U_star*/
/*    -Fill vector rhs*/
void computeRHS(data_Sim *sim, double *rhs, PetscInt rowStart, PetscInt rowEnd) {

    // YOU MUST IMPOSE A ZERO-MASS FLOW HERE ...
    // KESKE C KSA ???

    int i, j, start;
    int k = sim->ny + 2; // number of u values for each x position (nx + 1)
    int l = sim->ny + 1; // number of v values for each x position (nx + 2)

    int r = rowStart;
    
    for (i = 0; i < sim->nx; i++) {
        for (j = start; j < sim->ny; j++) {
            rhs[r++] = ((sim->u_star[(i+1)*k+j+1] - sim->u_star[i*k+j+1]) 
                      + (sim->v_star[(i+1)*l+j+1] - sim->v_star[(i+1)*l+j])) / sim->dt;
        }
    }

    rhs[rowStart] = 0.;  // set value of phi at inflow

        /*WRITE HERE (nabla dot u_star)/dt at each mesh point r*/
        /*Do not forget that the solution for the Poisson equation is defined within a constant.
        One point from Phi must then be set to an abritrary constant.*/

}

/*To call at each time step after computation of U_star. This function solves the poisson equation*/
/*and copies the solution of the equation into your vector Phi*/
/*More than probably, you should need to add arguments to the prototype ... */
/*Modification to do :*/
/*    - Change the call to computeRHS as you have to modify its prototype too*/
/*    - Copy solution of the equation into your vector PHI*/
void poisson_solver(data_Sim *sim, Poisson_data *data) {

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
    computeRHS(sim, rhs, rowStart, rowEnd); /*MODIFY THE PROTOTYPE HERE*/
    VecRestoreArray(b, &rhs);


    /*Solve the linear system of equations */
    KSPSolve(sles, b, x);
    KSPGetIterationNumber(sles, &its);
    PetscPrintf(PETSC_COMM_WORLD, "Solution to Poisson eqn in %d iterations \n", its);

    VecGetArray(x, &sol);

    int r;
    for(r=rowStart; r<rowEnd; r++){
        sim->phi[r] = sol[r];
        /*YOUR VECTOR PHI[...]*/ // = sol[r];
    }

    VecRestoreArray(x, &sol);
    free(sol);
    free(rhs);
}

/*This function is called only once during the simulation, i.e. in initialize_poisson_solver.*/
/*In its current state, it inserts unity on the main diagonal.*/
/*More than probably, you should need to add arguments to the prototype ... .*/
/*Modification to do in this function : */
/*   -Insert the correct factor in matrix A*/
void computeLaplacianMatrix(data_Sim *sim, Mat A, int rowStart, int rowEnd) {
    // ugly "if" in the loops, but only called once, so who cares 

    int i, j, idx, start, mask;
    int k = sim->ny;
    double count;
    
    MatSetValue(A, 0, 0, 1.0, INSERT_VALUES); // set reference pressure in origin box (i = 0, j = 0)

    for (i = 0; i < sim->nx; i++) {
        
        start = (i == 0) ? 1 : 0; // start at 1 if (i == 0), start at 0 otherwise

        for (j = start; j < sim->ny; j++) {
            idx = i * k + j;
            count = 0.;

            mask = (D_IN * sim->n <= i) && (i < (D_IN + LBOX) * sim->n) && (D_BOT * sim->n <= j) && (j < (D_BOT + 1) * sim->n);
            if (mask) { // inside rectangle
                count = 1;  // set diagonal to 1
            }
            else {
                if ((i != 0) && (i != (D_IN + LBOX) * sim->n)) { // it has left neighbor
                    MatSetValue(A, idx, idx-k, 1.0, INSERT_VALUES);
                    count += 1.;
                }
                if ((i != sim->nx - 1) && (i != D_IN * sim->n - 1)) { // it has right neighbor
                    MatSetValue(A, idx, idx+k, 1.0, INSERT_VALUES);
                    count += 1.;
                }
                if ((j != 0) && (j != (D_BOT + 1) * sim->n)) { // it has neighbor below
                    MatSetValue(A, idx, idx-1, 1.0, INSERT_VALUES);
                    count += 1.;
                }
                if ((j != sim->ny - 1) && (j != D_BOT * sim->n - 1)) { // it has neighbor above
                    MatSetValue(A, idx, idx+1, 1.0, INSERT_VALUES);
                    count += 1.;
                }
            }
            
            MatSetValue(A, idx, idx, -count, INSERT_VALUES);  // diagonal
        }
    }
    
    /*USING MATSETVALUE FUNCTION, INSERT THE GOOD FACTOR AT THE GOOD PLACE*/
    /*Be careful; the solution from the system solved is defined within a constant.
    One point from Phi must then be set to an abritrary constant.*/

}

/*To call during the initialization of your solver, before the begin of the time loop*/
/*Maybe you should need to add an argument to specify the number of unknows*/
/*Modification to do in this function :*/
/*   -Specify the number of unknows*/
/*   -Specify the number of non-zero diagonals in the sparse matrix*/
PetscErrorCode initialize_poisson_solver(data_Sim *sim, Poisson_data *data) {
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

    computeLaplacianMatrix(sim, data->A, rowStart, rowEnd);
    ierr = MatAssemblyBegin(data->A, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);
    ierr = MatAssemblyEnd(data->A, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);

    /* Create the Krylov context */
    KSPCreate(PETSC_COMM_WORLD, &(data->sles));
    KSPSetOperators(data->sles, data->A, data->A);
    KSPSetType(data->sles, KSPGMRES); // KSPGMRES seems the best, it will not be used if PC LU.
    PC prec;
    KSPGetPC(data->sles, &prec);
    PCSetType(prec, PCLU);
    // KSPSetFromOptions(data->sles); // to uncomment if we want to specify the solver to use in command line. Ex: mpirun -ksp_type gmres -pc_type gamg
    KSPSetTolerances(data->sles, 1.e-12, 1e-12, PETSC_DEFAULT, PETSC_DEFAULT);
    KSPSetReusePreconditioner(data->sles, PETSC_TRUE);
    KSPSetUseFischerGuess(data->sles, 1, 4);
    KSPGMRESSetPreAllocateVectors(data->sles);

    PetscPrintf(PETSC_COMM_WORLD, "Assembly of Mattrix and Vectors is done \n");

    return ierr;
}

/*To call after the simulation to free the vectors needed for Poisson equation*/
/*Modification to do : nothing */
void free_poisson_solver(Poisson_data* data) {
    MatDestroy(&(data->A));
    VecDestroy(&(data->b));
    VecDestroy(&(data->x));
    KSPDestroy(&(data->sles));
}
