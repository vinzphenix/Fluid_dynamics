//
//  thomas.h
//  MECA2660_hw 2016
//
//  Created by Thomas Gillis on 16/02/16.
//  Copyright © 2016 Thomas Gillis. All rights reserved.
//

#ifndef thomas_h
#define thomas_h

#include <stdlib.h>

void solve_Ac_thomas(const int n, const double a, const double b, const double c,
                     double *x, double *q, double *a1) {
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
    a1[0] = c / a;
    q[0] = q[0] / a;
    for (i = 1; i < n; ++i) {
        a1[i] = c / (a - b * a1[i - 1]);
        q[i] = (q[i] - b * q[i - 1]) / (a - b * a1[i - 1]);
    }

    // backward pass
    x[n - 1] = q[n - 1];
    for (i = n - 2; i >= 0; --i) {
        x[i] = q[i] - a1[i] * x[i + 1];
    }
}

void solve_period_3diag(const int n, const double a, const double b, const double c,
                        double *x, double *q) {
    /*
     This function solve a constant tridiagonal problem using Thomas algorithm.
     The matrix has the following structure:
     -                 -
     |a c             b|
     |b a c            |
     |  b a c          | = A        such that A x = q
     |        ...      |
     |c             b a|
     -                 -
     n: is the size of the problem
     a: is the main diagonal
     b: is the lower diagonal
     c: is the upper diagonal
     q: the rhs /!\ it overwrites the array /!\
     x: the solution is returned in this array /!\ it overwrites the array /!\
     */
    int i;

    double *q2 = (double *)calloc(n - 1, sizeof(double));  // needs to have zeros !
    double *x1 = (double *)malloc((n - 1) * sizeof(double));
    double *x2 = (double *)malloc((n - 1) * sizeof(double));
    double *a1 = (double *)malloc((n - 1) * sizeof(double));

    // boundary points treatment
    q2[0] = -b;
    q2[n - 2] = -c;

    // solve the different x1 and x2
    solve_Ac_thomas(n - 1, a, b, c, x1, q, a1);
    solve_Ac_thomas(n - 1, a, b, c, x2, q2, a1);

    // compute x_n
    x[n - 1] = (q[n - 1] - c * x1[0] - b * x1[n - 2]) / (a + c * x2[0] + b * x2[n - 2]);

    // compute the solution
    for (i = 0; i < n - 1; ++i) {
        x[i] = x1[i] + x2[i] * x[n - 1];
    }

    free(q2);
    free(x1);
    free(x2);
    free(a1);
}

#endif /* thomas_h */