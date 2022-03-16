#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>

#include "util.h"

int main(){

    /* parameters */

    double L = 20.0;
    double dx = 1e-3;
    
    int dim = (int)(2.0*L/dx);

    ParamF par;
    par.E = 0.8;
    par.xi = 1.0;

    double k = sqrt(par.E / par.xi);

    double x[dim], v[dim];
    complex double psi[dim];

    /* initial conditions */
    psi[0] = cexp(I*k*(-L));
    psi[1] = cexp(I*k*(-L+dx));

    /* solve the equation */
    fill_potential(v,x,dim);
    solve_numerov(x,psi,L,dx,F_step,&par);

    /* save the data */
    FILE *file;
    file = fopen("solution.csv","w");
    save_data(x,v,psi,dim,file);
    fclose(file);



}