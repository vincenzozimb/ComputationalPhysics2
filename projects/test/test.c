#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>

#include "bisection.h"
#include "zeros_newton.h"
#include "print_routines.h"
#include "numerov.h"

/* typedef */
typedef struct Param{
    double E, L, h;
    int l;
}Param;

/* functions */
double F(double x, void *p);
void solve_numerov(double x[], double psi[], double h, double F(double, void*), int dim, void *p);
double shooting(double E, void *p);
void normalize(double psi[], double h, int dim);
double delta(double E, void *p);

/* ------------------------------------------- */
int main(){

    /* parameters */
    double L = 5.0;     // with this value go up to 10 bound states
    double h = 1e-3;

    int l = 0;

    double E = 0.1;
    double dE = 0.1;
    int n = 8;
    double eig[n];

    Param par = {E,L,h,l};



}
/* ------------------------------------------- */

/* functions */
double F(double x, void *p){

    Param* par = (Param *)p;
    double E = par->E;

    return x * x - 2.0 * E;

}

void solve_numerov(double x[], double psi[], double h, double F(double, void*), int dim, void *p){
    for(int i=2; i<dim; i++){
        x[i] = x[i-1] + h;
        psi[i] = numerov_step(x[i-1],h,psi[i-1],psi[i-2],F,p);
    }
}

double shooting(double E, void *p){

    /* extract parameters */
    Param *par = (Param*)p;
    par->E = E;
    double L = par->L;
    double h = par->h;

    /* define variables */
    int dim = (int)(2.0*L/h);
    double x[dim], psi[dim];

    /* initial conditions */
    x[0] = -L;
    x[1] = -L + h;
    psi[0] = exp(-x[0] * x[0] / 2);
    psi[1] = exp(-x[1] * x[1] / 2);

    /* evolve with Numerov algorithm */
    solve_numerov(x,psi,h,F,dim,p);
    normalize(psi,h,dim);

    // FILE *file;
    // file = fopen("data.csv","w");
    // fprint_two_vec(file,x,psi,dim);
    // fclose(file);


    /* return */
    return psi[dim-2];

}

void normalize(double psi[], double h, int dim){
    
    double norm = 0.0;
    for(int i=0; i<dim; i++){
        norm += psi[i] * psi[i];
    }
    norm *= h;
    norm = sqrt(norm);

    for(int i=0; i<dim; i++){
        psi[i] /= norm;
    }

}

double delta(double E, void *p){

    /* extract parameters */
    Param *par = (Param *)p;
    
    par->E = E;
    double L = par->L;
    double h = par->h;
    int l = par->l;

    /* calculation inversion point */
    double x0 = zero_newton(F,0.5,p);
    // printf("x0 = %lf\n",x0);

    /* define variables */
    int dimL = (int)(x0/h);
    int dimR = (int)((L-x0)/h);

    double xL[dimL], xR[dimR], psiL[dimL], psiR[dimR];

    /* initial conditions */
    xL[0] = h;
    xL[1] = xL[0] + h;

    xR[0] = L;
    xR[1] = xR[0] - h;

    psiL[0] = pow(xL[0],l+1);
    psiL[1] = pow(xL[1],l+1);

    psiR[0] = exp(-2.0 * E * xR[0]);
    psiR[1] = exp(-2.0 * E * xR[1]);


    /* evolve with Numerov */
    


}