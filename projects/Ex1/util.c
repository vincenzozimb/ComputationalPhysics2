#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_sf_bessel.h>

#include "util.h"
#include "numerov.h"

/* functions */
double lj(double r){
    return 4 * ( pow(r,-12) - pow(r,-6) );
}

double F(double r, void *par){
    Param *p = (Param *)par;
    double xi = p->xi;
    double E = p->E;
    int l = p->l;
    
    return (1/xi) * ( lj(r) - E ) + (double)l*(l+1)/(r*r);
}

void solve_numerov(double r[], double u[], int dim, double dr, double F(double, void*), void *p){
    int i = 2;
    while(i < dim){
        u[i] = numerov_step(r[i-1],dr,u[i-1],u[i-2],F,p);
        r[i] = r[i-1] + dr;
        i++;
    }
}

double phase_shift(double k, int l, double r1, double r2, double u1, double u2){
    double K = (r1 * u2) / (r2 * u1);
    return atan( (gsl_sf_bessel_jl(l,k*r2) - K*gsl_sf_bessel_jl(l,k*r1)) / ( gsl_sf_bessel_yl(l,k*r2) - K*gsl_sf_bessel_yl(l,k*r1) ) );
}


