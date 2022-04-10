#include <stdio.h>
#include <stdlib.h>

int stupid() {

    int nx = 5;
    int ny = 3;

    double *u_data = calloc((nx+1) * (ny+2), sizeof(double));
    double **u = malloc((nx + 1) * sizeof(double *));


    printf("starting\n");
    for (int i = 0 ; i <= nx; i++) {
        printf("creation: row %d\n", i);
        u[i] = u_data + i * (ny+2);
    }

    u[1][1] = 7.;
    u[2][2] = 8.;
    u[3][0] = 9.;

    // corners
    u[0][0] = 1.;
    u[0][4] = 2.;
    u[5][0] = 3.;
    u[5][4] = 4.;

    for (int j = ny+1; j >= 0; j--) {
        for (int i = 0 ; i < nx+1; i++) {
            printf("%5.2lf\t", u[i][j]);
        }
        printf("\n");
    }

    free(u);
    free(u_data);

    return 0;
}
