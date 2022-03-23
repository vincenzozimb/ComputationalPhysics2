#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

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

double Fho(double r, void *par){
    Param *p = (Param *)par;
    double E = p->E;
    
    return r * r- 2.0 * E;
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

double delta(double E, void *p){
Param *par = (Param *)p;
    double L = par->L;
    double h = par->h;
    double x0 = sqrt(2.0 * E);
    int n = par->n;
    assert(n >= 0); // necessary to debug

    int dimL = (int)(L+x0)/h;
    int dimR = (int)(L-x0)/h;

    double xL[dimL], xR[dimR];
    double phiL[dimL], phiR[dimR];

    par->E = E;
    par->x0 = x0;

    /* ----- Initial conditions and evolution ----- */
    xL[0] = - L/2.0;     xL[1] = xL[0] + h;
    xR[0] = L/2.0;     xR[1] = xR[0] - h;
    phiL[0] = pow(-xL[0], n) * exp(- xL[0]*xL[0] /2.0);
    phiL[1] = pow(-xL[1], n) * exp(- xL[1]*xL[1] /2.0);
    initialize_r(xL, h, dimL);
    solve_numerov(xL, phiL, dimL, h, Fho, par);


    phiR[0] = pow(-xR[0], n) * exp(- xR[0]*xR[0] /2.0);
    phiR[1] = pow(-xR[1], n) * exp(- xR[1]*xR[1] /2.0);
    initialize_r(xR, -h, dimR);
    solve_numerov(xR, phiR, dimR, -h, Fho, par);

    /* ----- Adjust continuity ----- */
    double R = phiL[dimL - 1]/phiR[dimR -1];
    for (int ii = 0; ii < dimR; ii++)
    {
        phiR[ii] *= R;
    }
    
    /* ----- Return value ----- */

    double ris = phiL[dimL-2] + phiR[dimR-2];
    ris -= (2.0 + h * h * Fho(x0, par)) * phiL[dimL-1];
    ris /= h;

    return ris;
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