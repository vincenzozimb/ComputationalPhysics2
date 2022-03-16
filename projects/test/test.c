#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>

#include "numerov.h"
#include "print_routines.h"

/* typedef */
typedef struct ParamF{
    double xi,E;
}ParamF;

/* functions */
double potential(double x);
double F(double x, void *p);
void solve_numerov(double x[], double v[], complex double psi[], double L, double dx, double F(double, void *), void *p, FILE *outfile);



/* ------------------------------------------- */
int main(){

    /* parameter */
    double L = 10.0;
    double dx = 1e-3;

    int dim = (int)(2.0 * L / dx);

    double x[dim], v[dim];
    complex double psi[dim];

    ParamF par;
    par.xi = 1.0;
    par.E = 1.2;

    double k = sqrt(par.E / par.xi);

    /* initial conditions */
    psi[0] = cexp(-I*k*(-L));
    psi[1] = cexp(-I*k*(-L+dx));

    /* solve the equation */
    FILE *file;
    file = fopen("solution.csv","w");
    solve_numerov(x,v,psi,L,dx,F,&par,file);
    fclose(file);

}
/* ------------------------------------------- */

/* functions */
double potential(double x){
    if(fabs(x) < 0.5){
        return 1.0;
    }else{
        return 0.0;
    }
}

double F(double x, void *p){

    ParamF *par = (ParamF *)p;
    double xi = par->xi;
    double E = par->E;

    return (potential(x) - E) / xi;

}

void solve_numerov(double x[], double v[], complex double psi[], double L, double dx, double F(double, void *), void *p, FILE *outfile){

    /* assuming the first two values of psi to be already initialized */
    int dim = (int)(2.0 * L / dx);

    x[0] = -L;
    x[1] = x[0] + dx;
    v[0] = potential(x[0]);
    v[1] = potential(x[1]);

    for(int i=2; i<dim; i++){
        x[i] = x[i-1] + dx;
        v[i] = potential(x[i]);
        psi[i] = numerov_step(x[i-1],dx,psi[i-1],psi[i-2],F,p);
    }

    /* saving data */
    for(int i=0; i<dim; i++){
        fprint_double(outfile,x[i]);
        fprint_double(outfile,v[i]);
        fprint_double_newline(outfile,creal(psi[i]));
    }


}

