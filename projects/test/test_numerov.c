#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>

#include "bisection.h"
#include "zeros_newton.h"
#include "print_routines.h"
#include "numerov.h"

/* typedef */
typedef struct ParamPot{
    double R, rho, E;
    int l;
}ParamPot;

/* functions */
double F(double r, void *p);
double F_ho(double r, void *p);
void solve_numerov(double x[], double psi[], double h, double F(double, void*), int dim, void *p);
// double shooting(double E, void *p);
void normalize(double psi[], double h, int dim);
void print_wf(double r[], double psi[], int dim);

/* ------------------------------------------- */
int main(){

    /* system parameters */
    int N = 20; // number of electrons
    double rs = 3.93; // 3.93 for Na and 4.86 for K

    double R = rs * pow((double)N,1.0/3.0); // radius of the cluster (of its harmonic part)
    double rho = 3.0 / (4.0 * M_PI * rs * rs * rs); // density of the jellium

    /* code parameters */
    double L = 3.0 * R; // infrared cutoff. Space interval goes from 0 from L
    double h = 1e-3;; // ultraviolet cutoff
    int dim = (int)(L/h); // number of discretizations in the finite difference method

    double E = 1.5;
    int l = 0;

    double r[dim], psi[dim];
    ParamPot par = {R,rho,E,l};

    r[0] = h;
    r[1] = r[0] + h;
    psi[0] = pow(r[0],l+1);
    psi[1] = pow(r[1],l+1);

    solve_numerov(r,psi,h,F_ho,dim,&par);
    normalize(psi,h,dim);
    print_wf(r,psi,dim);


}
/* ------------------------------------------- */

/* functions */
double potential(double r, void *p){

    ParamPot *par = (ParamPot *)p;

    double R = par->R;
    double rho = par->rho;

    double pot;
    if(r <= R){
        pot = r * r / 3.0 - R * R;
    }else{
        pot = -2.0/3.0 * R * R * R / r; 
    }

    pot *= 2.0 * M_PI * rho;
    return pot;

}

double F(double r, void *p){

    ParamPot* par = (ParamPot *)p;
    int l = par->l;
    double E = par->E;

    return 2.0 * (potential(r,p) + (double)l*(l+1)/(2.0*r*r) - E);

}

double F_ho(double r, void *p){

    ParamPot* par = (ParamPot *)p;
    int l = par->l;
    double E = par->E;

    return 2.0 * (0.5*r*r + (double)l*(l+1)/(2.0*r*r) - E);

}

void solve_numerov(double x[], double psi[], double h, double F(double, void*), int dim, void *p){
    for(int i=2; i<dim; i++){
        x[i] = x[i-1] + h;
        psi[i] = numerov_step(x[i-1],h,psi[i-1],psi[i-2],F,p);
    }
}

// double shooting(double E, void *p){

//     /* extract parameters */
//     Param *par = (Param*)p;
//     par->E = E;
//     double L = par->L;
//     double h = par->h;

//     /* define variables */
//     int dim = (int)(2.0*L/h);
//     double x[dim], psi[dim];

//     /* initial conditions */
//     x[0] = -L;
//     x[1] = -L + h;
//     psi[0] = exp(-x[0] * x[0] / 2);
//     psi[1] = exp(-x[1] * x[1] / 2);

//     /* evolve with Numerov algorithm */
//     solve_numerov(x,psi,h,F,dim,p);
//     normalize(psi,h,dim);

//     // FILE *file;
//     // file = fopen("data.csv","w");
//     // fprint_two_vec(file,x,psi,dim);
//     // fclose(file);


//     /* return */
//     return psi[dim-2];

// }

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

void print_wf(double r[], double psi[], int dim){

    FILE *file;
    file = fopen("wf.csv","w");
    fprint_two_vec(file,r,psi,dim);
    fclose(file);

}