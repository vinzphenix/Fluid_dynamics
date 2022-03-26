#ifndef _PROJECT_H_
#define _PROJECT_H_

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define RE 500.
#define CFL 1.
#define FOURIER 1.


/*
 x_ = x / Hbox
 y_ = y / Hbox
 u_ = u / U0
 v_ = v / U0
 t_ = t / (Hbox/U0)
 p_ = p / (rho U0^2)
*/

typedef struct {
    int nx, ny;
    double h, dt;
    double *u;
    double *v;
    double *p;
} data_Sim;

#endif

