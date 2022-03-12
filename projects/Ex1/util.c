#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

#include "util.h"
#include "numerov.h"
#include "print_routines.h"

/* functions */
double potential(double x){
    if(fabs(x) <= 1.0){
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

void initialize_pot(double x[], double V[], int dim){
    for(int i=0; i<dim; i++){
        V[i] = potential(x[i]);
    }
}

void save_data(double x[], double V[], complex double psi[], int dim, FILE *outfile){
    for(int i=0; i<dim; i++){
        fprint_double(outfile,x[i]);
        fprint_double(outfile,V[i]);
        fprint_double(outfile,creal(psi[i]));
        fprint_double_newline(outfile,cimag(psi[i]));
    }
}



