#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

#include "util.h"
#include "numerov.h"
#include "print_routines.h"

/* functions */
double V(double x){
    //return exp(-x * x / 2.0);
    if (fabs(x) < 0.5) {
        return 1.0;
    }
    else return 0.0;
}

double F(double x, void *p){
    
    Param *par = (Param *)p;
    double xi = par->xi;
    double E = par->E;

    return (1/xi) * (V(x) - E);
}

void solve_numerov(double x[], complex double psi[], int dim, double dx, double F (double, void *p), void *p, FILE *outfile){

    /* assuming that the first two values of x[] and psi[] are initialized */

    double pot[dim];
    pot[0] = V(x[0]);
    pot[1] = V(x[1]);

    int i = 2;
    while(i < dim){
        /* evolve */
        psi[i] = numerov_step(x[i-1],dx,psi[i-1],psi[i-2],F,p);
        x[i] = x[i-1] + dx;
        pot[i] = V(x[i]);      
        i++;
    }
    for(i=0;i<dim;i++){
        /* print */
        fprint_double(outfile, x[i]);
        fprint_double(outfile,pot[i]);
        fprint_double(outfile, creal(psi[i]));
        fprint_double(outfile, cimag(psi[i]));
        fprint_double_newline(outfile, (double)(psi[i] * conj(psi[i])));
    }
}

double normalize(complex double psi[], int dim, double dx){
    double norm = 0.0;
    for(int i=0; i<dim; i++){
        norm += psi[i] * conj(psi[i]) * dx;
    }
    for(int i=0; i<dim; i++){
        psi[i] /= norm;
    }
    return norm;
}


