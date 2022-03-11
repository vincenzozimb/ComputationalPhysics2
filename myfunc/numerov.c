#include "numerov.h"
#include <complex.h>

double complex numerov_step(double x, double dx,
                            complex double psi_curr, complex double psi_prec,
                            double F(double, void *), void *p){

    complex double ris;
    ris = (2.0 + 5.0 * dx * dx * F(x, p) / 6.0) * psi_curr;
    ris -= (1.0 - dx * dx * F(x - dx, p) / 12.0) * psi_prec;
    ris /= 1.0 - dx * dx * F(x + dx, p) / 12.0;
    return ris;
}

void execute_numerov_complex(double x[], double dx, complex double psi[], int dim,
                             double F(double, void *), void *p){
    int i = 2;
    while (i < dim){
        psi[i] = numerov_step(x[i-1], dx, psi[i-1], psi[i-2], F, p);
        x[i] = x[i-1] + dx;
        i++;
    }
}