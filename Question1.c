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


int main() {

}