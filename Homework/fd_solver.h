#ifndef fd_solver_h
#define fd_solver_h

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define SAVE 1
#define PRINT(x) x

//char *filename = "./data/solution_E2_64.txt";
//char *filename = "./data/solution.txt";
char *path = "./data/solution";
char filename[50];

// Simulation parameters
#define N 64
#define TEND 1.
#define L 1.
#define C 1.
#define SIGMA (1. / 16.)
#define UMAX 1.

#define SCHEME_A 'E'
#define SCHEME_B '4'

#define CFL 1.4
// CFL = (E2: 2.828) (E4: 2.061) (E6: 1.783) (I4: 1.632) (I6: 1.421)


double ALPHA[4] = {0, 0.5, 0.5, 1.};
double GAMMA[4] = {0.16666666666, 0.33333333333, 0.33333333333, 0.16666666666};

/*typedef struct data_Sim_alias {
    int N, M;
    double c, sigma, U_max, L, h, dt, CFL;
    double *u, *us, *ul, *du, *q;
    char *scheme;
    void (*f_eval)(struct data_Sim_alias *);
} data_Sim;*/

typedef struct data_Sim_alias {
    int M;
    double h, dt;
    double *u, *us, *ul, *du, *q;
} data_Sim;

void f_eval(data_Sim *sim);

#endif