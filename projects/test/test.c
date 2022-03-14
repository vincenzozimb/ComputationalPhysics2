#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>

#include "util.h"
#include "print_routines.h"

int main(){

    /* parameters */

    int dim = 1000;
    
    double x[dim], v[dim];
    complex double psi[dim];
    double L = 10.0;
    double dx = L / (dim-1);

    ParamF par;
    par.xi = 1.0;
    par.E = 0.8;

    double k = sqrt(par.E / par.xi);

    /* initial conditions */
    x[0] = -L / 2.0;
    x[1] = x[0] + dx;
    psi[0] = cexp(I*k*x[0]);
    psi[1] = cexp(I*k*x[1]);

    /* solving the equation */
    solve_numerov(x,v,psi,dim,dx,F,&par);
    
    /* saving the data */
    FILE *file;
    file = fopen("solution.csv","w");
    save_data(x,v,psi,dim,file); 
    fclose(file);



}