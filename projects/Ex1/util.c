#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "util.h"
#include "numerov.h"
#include "bessel_func.h"

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

void initialize_r(double r[], double dr, int dim){
    /* assuming the first two values of r[] to be already initialized */
    for(int i=2; i<dim; i++){
        r[i] = r[i-1] + dr;
    }
}

void solve_numerov(double r[], double u[], int dim, double dr, double F(double, void*), void *p){
    /* assuming the first two values of u[] to be already initialized */
    for(int i=2; i<dim; i++){
        u[i] = numerov_step(r[i],dr,u[i-1],u[i-2],F,p);
    }
}

double phase_shift(double r1, double r2, double u1, double u2, double k, int l, FuncBessel *f1, FuncBessel *f2){

    double K = (r1 * u2) / (r2 * u1);
    double delta = atan( (f2->j_curr - K * f1->j_curr) / (f2->n_curr - K * f1->n_curr) );


    /* uptade the Bessel functions */
    double j1 = f1->j_prec;
    double n1 = f1->n_prec;
    double j2 = f2->j_prec;
    double n2 = f2->n_prec;

    f1->j_prec = f1->j_curr;
    f1->n_prec = f1->n_curr;
    f2->j_prec = f2->j_curr;
    f2->n_prec = f2->n_curr;

    f1->j_curr = recursive_bessel(k*r1,l,j1,f1->j_prec);
    f1->n_curr = recursive_bessel(k*r1,l,n1,f1->n_prec);
    f2->j_curr = recursive_bessel(k*r2,l,j2,f2->j_prec);
    f2->n_curr = recursive_bessel(k*r2,l,n2,f2->n_prec);
    

    /* return */
    return delta;

}