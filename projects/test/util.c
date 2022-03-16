#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

#include "util.h"
#include "numerov.h"
#include "print_routines.h"

/* functions */
double potential(double x){
    
    if(x < 0.0){
        return 0.0;
    }else{
        return 1.0;
    }

}

double F_step(double x, void *p){

    ParamF *par = (ParamF *)p;
    double E = par->E;
    double xi = par->xi;

    double pot;
    if(x < 0.0){
        pot = 0.0;
    }else{
        pot = 1.0;
    }

    return (pot - E) / xi;

}

void solve_numerov(double x[], complex double psi[], double L, double dx, double F(double, void *), void *p){

    int dim = (int)(2.0*L/dx);

    x[0] = -L;
    x[1] = -L + dx;

    /* assuming the first two values of psi to be already initialized */
    for(int i=2; i<dim; i++){
        psi[i] = numerov_step(x[i-1],dx,psi[i-1],psi[i-1],F,p);
        x[i] = x[i-1] + dx;
    }

}

void fill_potential(double v[], double x[], int dim){
    
    for(int i=0; i<dim; i++){
        v[i] = potential(x[i]);
    }

}

void save_data(double x[], double v[], complex double psi[], int dim, FILE *outfile){

    for(int i=0; i<dim; i++){
        fprint_double(outfile,x[i]);
        fprint_double(outfile,v[i]);
        fprint_double_newline(outfile,creal(psi[i]));
    }

}
