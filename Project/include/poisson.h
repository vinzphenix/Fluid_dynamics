#ifndef _POISSON_H_
#define _POISSON_H_

/*To include in the file in which you will call initialize_poisson_solver and poisson_solver*/

#include <petsc.h>
#include <petscsys.h>
#include "project.h"
#include "mpi.h"

//Structure storing petsc vectors

typedef struct {

	Vec b;
	Vec x;
	Mat A;
	KSP sles;

} Poisson_data;

PetscErrorCode initialize_poisson_solver(Sim_data *sim, Poisson_data* data);
int poisson_solver(Sim_data *sim, Poisson_data *data);
void free_poisson_solver(Poisson_data* data);
void test_poisson(Sim_data *sim, Poisson_data *poisson);

#endif
