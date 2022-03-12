#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

#include "util.h"
#include "numerov.h"
#include "print_routines.h"

/* functions */
double F_gauss(double x, void *p){
    
    Param *par = (Param *)p;
    double xi = par->xi;
    double E = par->E;

    return (1/xi) * (-exp(-x*x/2) - E);
}

void solve_numerov(double x[], complex double psi[], int dim, double dx, double F (double, void *p), void *p, FILE *outfile){

    /* assuming that the first two values of x[] and psi[] are initialized */

    complex double psi_new;
    int i = 2;
    while(i < dim){
        /* evolve */
        psi_new = numerov_step(x[i-1],dx,psi[i-1],psi[i-2],F,p);
        psi[i] = psi_new;
        x[i] = x[i-1] + dx;
        /* print */
        fprint_double(outfile, x[i]);
        fprint_double(outfile, creal(psi[i]));
        fprint_double(outfile, cimag(psi[i]));
        fprintf(outfile, "\n");
        
        i++;
    }


}



