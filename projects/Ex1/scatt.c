#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <complex.h>
#include <math.h>

#include "util.h"
#include "numerov.h"
#include "print_routines.h"

int main(){

    /* first I will try to solve the SE with numerov to remember how to do it */

    /* problem parameters */
    Param par;
    par.xi = 1.0;
    par.E = 1.2;
    
    double k = sqrt(par.E / par.xi);
    double L = 2 * 10.0;
    double dx = 0.001;

    int dim = (int)(L / dx);
    double x[dim];
    complex double psi[dim];

    /* initial conditions */
    x[0] = -L / 2.0;
    x[1] = x[0] + dx;
    psi[0] = cexp(I*k*(x[0]));
    psi[1] = cexp(I*k*(x[1]));

    /* solve the equation */
    FILE *file;
    file = fopen("solution.csv", "w");

    solve_numerov(x,psi,dim,dx,F,&par,file);

    fclose(file);



}