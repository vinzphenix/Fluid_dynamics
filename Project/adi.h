//
//  thomas.h
//  MECA2660_hw 2016
//
//  Created by Thomas Gillis on 16/02/16.
//  Copyright © 2016 Thomas Gillis. All rights reserved.
//

#ifndef thomas_h
#define thomas_h

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

typedef struct {
    int nx, ny, size;
    double *a, *b, *c;
    double *q;
} ADI_data;


void solve_Ac_thomas(const int n, const double a, const double b, const double c,
                     double *x, double *q, double *at) {
    /*
     This function solve a constant tridiagonal problem using Thomas algorithm.
     The matrix has the following structure:
     -                 -
     |a c              |
     |b a c            |
     |  b a c          | = Ac       such that Ac x = q
     |        ...      |
     |            b a  |
     -                 -
     n: is the size of the problem
     a: is the main diagonal
     b: is the lower diagonal
     c: is the upper diagonal
     q: the rhs /!\ it overwrites the array /!\
     x: the solution is returned in this array /!\ it overwrites the array /!\
     */

    int i;

    // forward pass
    at[0] = c / a;
    q[0] = q[0] / a;
    for (i = 1; i < n; ++i) {
        at[i] = c / (a - b * at[i - 1]);
        q[i] = (q[i] - b * q[i - 1]) / (a - b * at[i - 1]);
    }

    // backward pass
    x[n - 1] = q[n - 1];
    for (i = n - 2; i >= 0; --i) {
        x[i] = q[i] - at[i] * x[i + 1];
    }
}

// void solve_period_3diag(const int n, const double a, const double b, const double c,
//                         double *x, double *q, double *x1, double *at) {
//     /*
//      This function solve a constant tridiagonal problem using Thomas algorithm.
//      The matrix has the following structure:
//      -                 -
//      |a c             b|
//      |b a c            |
//      |  b a c          | = A        such that A x = q
//      |        ...      |
//      |c             b a|
//      -                 -
//      n: is the size of the problem
//      a: is the main diagonal
//      b: is the lower diagonal
//      c: is the upper diagonal
//      q: the rhs /!\ it overwrites the array /!\
//      x: the solution is returned in this array /!\ it overwrites the array /!\
//      */
//     int i;

//     // solve for the first x
//     solve_Ac_thomas(n - 1, a, b, c, x1, q, at);

//     // solve for the second x: reset and setup q
//     memset(q, 0, (n - 1) * sizeof(double));
//     q[0] = -b;
//     q[n - 2] = -c;
//     solve_Ac_thomas(n - 1, a, b, c, x, q, at);

//     // compute x_n
//     x[n - 1] = (q[n - 1] - c * x1[0] - b * x1[n - 2]) / (a + c * x[0] + b * x[n - 2]);

//     // compute the solution
//     for (i = 0; i < n - 1; ++i) {
//         x[i] = x1[i] + x[i] * x[n - 1];
//     }
// }

#endif /* thomas_h */
