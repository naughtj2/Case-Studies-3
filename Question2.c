#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define IDX(i, j, N) ((i) * (N) + (j))
#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

double f(double x1, double x2) {
    return 2 * M_PI * M_PI * sin(M_PI * x1) * sin(M_PI * x2);
}

// Setup for the poisson problem
void poisson_setup(int N, int size, double **A, double *b, double h) {
    *A = (double *) calloc(size * size, sizeof(double));
    if (*A == NULL) {
        perror("Matrix allocation failed");
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < N; ++i) {
        double x1 = (i + 1) * h;
        for (int j = 0; j < N; ++j) {
            double x2 = (j + 1) * h;
            int row = IDX(i, j, N);
            
            b[row] = f(x1, x2);

            (*A)[row * size + row] = 4.0;
            if (i > 0)
                (*A)[row * size + IDX(i - 1, j, N)] = -1.0;
            if (j < N - 1)
                (*A)[row * size + IDX(i + 1, j, N)] = -1.0;
            if (j > 0)
                (*A)[row * size + IDX(i, j - 1, N)] = -1.0;
            if (j < N - 1)
                (*A)[row * size + IDX(i, j + 1, N)] = -1.0;
            
            b[row] *= h * h;
        }

    }
}

// Apply A*y 
void applyA(int N, const double *y, double *Ay) {
    double h = 1.0 / (N + 1);
    double h2_inv = 1.0 / (h * h);

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            int k = IDX(i, j, N);
            
            Ay[k] = 4.0 * y[k];
            if (i > 0)
                Ay[k] -= y[IDX(i - 1, j, N)];
            if (i < N - 1)
                Ay[k] -= y[IDX(i + 1, j, N)];
            if (j > 0)
                Ay[k] -= y[IDX(i, j - 1, N)];
            if (j < N - 1)
                Ay[k] -= y[IDX(i, j + 1, N)];

            Ay[k] *= h2_inv;
        }
    }
}

// Dot Product
double dot(const double *a, const double *b, int n) {
    double sum = 0.0;
    for (int i = 0; i < n; ++i) 
        sum += a[i] * b[i];
    return sum;
}

