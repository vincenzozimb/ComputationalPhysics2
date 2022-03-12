#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <complex.h>
#include <math.h>

#include "util.h"
#include "numerov.h"

int main(){

    /* first I will try to solve the SE with numerov to remember how to do it */

    /* problem parameters */
    Param par;
    par.xi = 1.0;
    par.E = 1.2;
    
    double k = sqrt(par.xi)*sqrt(par.E);
    double L = 2 * 5.0;
    double dx = 0.01;

    int dim = (int)(L / dx);
    double x[dim];
    complex double psi[dim];

    /* initial conditions */
    x[0] = -L;
    x[1] = x[0] + dx;
    psi[0] = cexp(-I*k*(-L));
    psi[1] = cexp(-I*k*(-L + dx));

    /* solve the equation */
    FILE *file;
    file = fopen("solution.csv", "w");

    solve_numerov(x,psi,dim,dx,F_gauss,&par,file);



}