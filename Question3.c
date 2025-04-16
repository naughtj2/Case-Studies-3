#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <time.h>

double reltol = 1e-8;

// Dot Product
double dot(const double *a, const double *b, int n) {
    double sum = 0.0;
    for (int i = 0; i < n; ++i) 
        sum += a[i] * b[i];
    return sum;
}

// Matric Vecotr Multiplication
void matvec(int N, const double *A, const double *x, double *y) {
    for (int i = 0; i < N; ++i) {
        y[i] = 0.0;
        for (int j = 0; j < N; ++j)
            y[i] += A[i * N + j] * x[j];
    }
}

// Conjuagte Gradient for Dense Linear System
int conjugate_gradient_dense(int N, double *x, const double *b, const double *A, double *residuals) {
    int max_iter = 10000;
    int k = 0;

    double *r = calloc(N, sizeof(double));
    double *p = calloc(N, sizeof(double));
    double *Ap = calloc(N, sizeof(double));

    matvec(N, A, x, Ap);
    for (int i = 0; i < N; ++i)
        r[i] = b[i] - Ap[i];

    double r0_norm = sqrt(dot(r, r, N));
    double tol = reltol * r0_norm;

    for (int i = 0; i < N; ++i)
        p[i] = r[i];

    double rs_old = dot(r, r, N);
    double rs_new;

    while (k < max_iter) {
        matvec(N, A, p, Ap);
        double alpha = rs_old / dot(p, Ap, N);

        for (int i = 0; i < N; ++i)
            x[i] += alpha * p[i];
        for (int i = 0; i < N; ++i)
            r[i] -= alpha * Ap[i];

        rs_new = dot(r, r, N);
        residuals[k] = sqrt(rs_new);
        if (residuals[k] <= tol)
            break;

        for (int i = 0; i < N; ++i)
            p[i] = r[i] + (rs_new / rs_old) * p[i];

        rs_old = rs_new;
        ++k;
    }

    free(r);
    free(p);
    free(Ap);
    return k;
}

int main() {
    reltol = sqrt(DBL_EPSILON);

    int Ns[] = {100, 1000, 10000};
    int num_Ns = sizeof(Ns) / sizeof(Ns[0]);

    for (int idx = 0; idx < num_Ns; ++idx) {
        int N = Ns[idx];
        printf("Running CG for N = %d...\n", N);

        double *A = malloc(N * N * sizeof(double));
        double *b = malloc(N * sizeof(double));
        double *x = calloc(N, sizeof(double));
        double *residuals = calloc(N, sizeof(double));
        
        for (int i = 0; i < N; ++i) {
            b[i] = 1.0;
            for (int j = 0; j < N; ++j)
                A[i * N + j] = (double)(N - abs(i - j)) / N;
        }

        clock_t start = clock();
        int iters = conjugate_gradient_dense(N, x, b, A, residuals);
        clock_t end = clock();
        double time_spent = (double)(end - start) / CLOCKS_PER_SEC;

        char fname[64];
        sprintf(fname, "residuals_%d.txt", N);
        FILE *fres = fopen(fname, "w");
        for (int k = 0; k < iters; ++k)
            fprintf(fres, "%d %.12e\n", k, residuals[k]);
        fclose(fres);

        printf("CG finished in %d iterations (%.3f sec)\n\n", iters, time_spent);

        free(A);
        free(b);
        free(x);
        free(residuals);
    }

    return 0;
}
