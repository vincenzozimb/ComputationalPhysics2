#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "util.h"
#include "numerov.h"
#include "bessel_func.h"
#include "print_routines.h"
#include "bisection.h"

/* functions */

// GENERIC
void initialize_position(double x[], double h, int dim){
    /* assuming the first two values of r[] to be already initialized */
    for (int i = 2; i < dim; i++){
        x[i] = x[i-1] + h;
    }
}

void solve_numerov(double x[], double psi[], int dim, double h, double F(double, void *), void *p){
    /* assuming the first two values of u[] to be already initialized */
    for (int i = 2; i < dim; i++){
        psi[i] = numerov_step(x[i-1], h, psi[i-1], psi[i-2], F, p);
    }
}

// SCATTERING 
double lj(double r){
    return 4 * (pow(r, -12) - pow(r, -6));
}

double F_lj(double r, void *p){
    Param_F *par = (Param_F *)p;
    double xi = par->xi;
    double E = par->E;
    int l = par->l;

    return (1 / xi) * (lj(r) - E) + (double)l * (l + 1) / (r * r);
}

double phase_shift(double r1, double r2, double u1, double u2, double k, int l, FuncBessel *f1, FuncBessel *f2){

    double K = (r1 * u2) / (r2 * u1);
    double delta = atan((f2->j_curr - K * f1->j_curr) / (f2->n_curr - K * f1->n_curr));

    /* uptade the Bessel functions */
    double j1 = f1->j_prec;
    double n1 = f1->n_prec;
    double j2 = f2->j_prec;
    double n2 = f2->n_prec;

    f1->j_prec = f1->j_curr;
    f1->n_prec = f1->n_curr;
    f2->j_prec = f2->j_curr;
    f2->n_prec = f2->n_curr;

    f1->j_curr = recursive_bessel(k * r1, l, j1, f1->j_prec);
    f1->n_curr = recursive_bessel(k * r1, l, n1, f1->n_prec);
    f2->j_curr = recursive_bessel(k * r2, l, j2, f2->j_prec);
    f2->n_curr = recursive_bessel(k * r2, l, n2, f2->n_prec);

    /* return */
    return delta;
}

// BOUND STATES
double F_ho(double x, void *p){
    Param_F *par = (Param_F *)p;
    double E = par->E;

    return x * x - 2.0 * E;
}

double F_ho3d(double x, void *p){
    Param_F *par = (Param_F *)p;
    double E = par->E;
    int l = par->l;

    return x * x + ( (double)l*(l+1) / (x * x) ) - 2.0 * E;
}

void normalize(double psi[], double dx, int dim){
    /* calculate norm */
    double norm = 0.0;
    
    for(int i=0; i<dim; i++){
        norm += psi[i] * psi[i];
    }
    norm *= dx;
    norm = sqrt(norm);

    /* normalize */
    for(int i=0; i<dim; i++){
        psi[i] /= norm;
    }
}

void initial_condition_ho(double x[], double psi[], double L, double h, int n){
    /* position */
    x[0] = -L;
    x[1] = x[0] + h;

    /* wave function */
    psi[0] = pow(x[0],n) * exp(-x[0] * x[0] / 2.0);
    psi[1] = pow(x[1],n) * exp(-x[1] * x[1] / 2.0);
}

void initial_condition_ho3d(double xL[], double xR[], double psiL[], double psiR[], double L, double h, int l){
    /* position */
    xL[0] = h;
    xL[1] = xL[0] + h;

    xR[0] = L;
    xR[1] = xR[0] - h;

    /* wave function */
    psiL[0] = pow(xL[0],l+1);
    psiL[1] = pow(xL[1],l+1);

    psiR[0] = exp(-xR[0] * xR[0] / 2.0);
    psiR[1] = exp(-xR[1] * xR[1] / 2.0);
}

double run_for_delta(double E, void *p){
    /* extract parameters */
    Param_delta *pd = (Param_delta *)p;
    
    double h = pd->h;
    double L = pd->L;
    int n = pd->n;

    double x0 = sqrt(2.0 * E) - h;
    pd->x0 = x0;

    /* define variables */
    Param_F pf;
    pf.E = E;

    int dimL = 1 + (int)ceil((L + x0) / h);
    int dimR = 1 + (int)floor((L - x0) / h);

    double xL[dimL], xR[dimR], psiL[dimL], psiR[dimR];

    /* initial conditions */
    initial_condition_ho(xL,psiL,L,h,n);
    initial_condition_ho(xR,psiR,-L,-h,n);

    /* evolution */
    initialize_position(xL,h,dimL);
    solve_numerov(xL,psiL,dimL,h,F_ho,&pf);

    initialize_position(xR,-h,dimR);
    solve_numerov(xR,psiR,dimR,-h,F_ho,&pf);

    /* normalize */
    normalize(psiL,h,dimL);
    normalize(psiR,h,dimR);

    /* adjust continuity */
    double R = psiL[dimL-1] / psiR[dimR-1];
    for(int i=0; i<dimR; i++){
        psiR[i] *= R;
    }

    /* return delta */
    assert( fabs(psiL[dimL-1] - psiR[dimR-1]) < EPS );
    double ris = psiL[dimL-2] + psiR[dimR-2];
    ris -= (2.0 + h * h * F_ho(x0,&pf)) * psiL[dimL-1];
    ris /= h;

    return ris;

}

double run_for_delta3d(double E, void *p){
    /* extract parameters */
    Param_delta *pd = (Param_delta *)p;
    
    double h = pd->h;
    double L = pd->L;
    int l = pd->l;
    double x0;
    double xmin = 0.0;

    if (l == 0)
    {
        xmin = 0.0;
    } else if (l == 1)
    {
        xmin = 1.2;
    } else if (l == 2){
        xmin = 1.6;
    } else if (l == 3){
        xmin = 1.9;
    } else if (l == 4){
        xmin = 2.2;
    }

    Param_F pf;
    pf.E = E; 
    pf.l = l;
    x0 = bisection(F_ho3d,xmin+h,L-h,&pf);
    // printf("E: %lf\t x0 : %lf\n", E, x0);
    pd->x0 = x0;

    /* define variables */
    int dimL, dimR;

    dimL = 1 + (int)floor(x0 / h);
    dimR = 1 + (int)floor((L - x0) / h);

    double xL[dimL], xR[dimR], psiL[dimL], psiR[dimR]; 

    /* initial conditions */
    initial_condition_ho3d(xL,xR,psiL,psiR,L,h,l); 

    /* evolution */
    initialize_position(xL,h,dimL);
    solve_numerov(xL,psiL,dimL,h,F_ho3d,&pf);

    initialize_position(xR,-h,dimR);
    solve_numerov(xR,psiR,dimR,-h,F_ho3d,&pf);

    /* normalize */
    normalize(psiL,h,dimL);
    normalize(psiR,h,dimR);

    /* adjust continuity */
    double R = psiL[dimL-1] / psiR[dimR-1];
    for(int i=0; i<dimR; i++){
        psiR[i] *= R;
    }

    /* return delta */
    assert( fabs(psiL[dimL-1] - psiR[dimR-1]) < EPS );
    double ris = psiL[dimL-2] + psiR[dimR-2];
    ris -= (2.0 + h * h * F_ho3d(x0,&pf)) * psiL[dimL-1];
    ris /= h;

    return ris;

}
