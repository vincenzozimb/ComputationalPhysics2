#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

#include "util.h"
#include "numerov.h"
#include "print_routines.h"

/* functions */
double pot_step(double x){
    if(x > 0.0){
        return 1.0;
    }
    else return 0.0;
}

double F(double v, void *p){

    ParamF *par = (ParamF *)p;
    double xi = par->xi;
    double E = par->E;

    return (v - E) / xi;
}

void solve_numerov(double x[], double v[], complex double psi[], int dim, double dx, double F(double, void *), void *p){

    /* assuming that the first two values of x[] and psi[] are already initialized */
    
    ParamF *par = (ParamF *)p;
    double xi = par->xi;
    double E = par->E;

    v[0] = xi * F(x[0],par) + E;
    v[1] = xi * F(x[1],par) + E;
    
    for(int i=2; i<dim; i++){
        psi[i] = numerov_step(x[i-1],dx,psi[i-1],psi[i-2],F,p);
        x[i] = x[i-1] + dx;
        v[i] = xi * F(x[i],par) + E;
    }
}

void save_data(double x[], double v[], complex double psi[], int dim, FILE *outfile){

    for(int i=0; i<dim; i++){
        fprint_double(outfile,x[i]);
        fprint_double(outfile,v[i]);
        fprint_double_newline(outfile,creal(psi[i]));
    }

}
