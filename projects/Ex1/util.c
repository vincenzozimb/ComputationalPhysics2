#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

#include "util.h"
#include "numerov.h"
#include "print_routines.h"

/* functions */
double potential(double x){
    if(fabs(x) <= 0.5){
        return 1.0;
    }
    else return 0.0;
}

double F(double x, void *p){
    Param *par = (Param *)p;
    double xi = par->xi;
    double E = par->E;

    return (1.0/xi) * (potential(x) - E);
}

void solve_numerov(double x[], complex double psi[], int dim, double dx, double F(double, void *p), void *p, FILE *outfile){

    double V[dim];
    V[0] = potential(x[0]);
    V[1] = potential(x[1]);

    for(int i=2; i<dim; i++){
        psi[i] = numerov_step(x[i-1],dx,psi[i-1],psi[i-2],F,p);
        x[i] = x[i-1] + dx;
        V[i] = potential(x[i]);
    }

    for(int i=0; i<dim; i++){
        fprint_double(outfile,x[i]);
        fprint_double(outfile,V[i]);
        fprint_double_newline(outfile,creal(psi[i]));
    }

}


