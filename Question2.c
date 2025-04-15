#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define IDX(i, j, N) ((i) * (N) + (j))
#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif
#define TOL 1e-8

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

// Conjugate Gradient Solver
int conjugate_gradient(int N, int size, double *x, const double *b) {
    int max_iter = 10000;
    int k = 0;

    double *r = calloc(size, sizeof(double));
    double *p = calloc(size, sizeof(double));
    double *Ap = calloc(size, sizeof(double));

    applyA(N, x, Ap);

    for (int i = 0; i < size; ++i)
        r[i] = b[i] - Ap[i];

    for (int i = 0; i < size; ++i)
        p[i] = r[i];

    double rs_old = dot(r, r, size);
    double rs_new;

    // printf("Initial residual norm: %.12f\n", sqrt(rs_old));

    do {
        applyA(N, p, Ap);
        double alpha = rs_old / dot(p, Ap, size);

        for (int i = 0; i < size; ++i)
            x[i] += alpha * p[i];
        for (int i = 0; i < size; ++i)
            r[i] -= alpha * Ap[i];

        rs_new = dot(r, r, size);
        ++k;

        if (sqrt(rs_new) < TOL)
            break;

        for (int i = 0; i < size; ++i)
            p[i] = r[i] + (rs_new / rs_old) * p[i];

        rs_old = rs_new;        
    } while (k < max_iter);

    free(r);
    free(p);
    free(Ap);

    // printf("Final residual norm: %.12f\n", sqrt(rs_new));
    return k;
}