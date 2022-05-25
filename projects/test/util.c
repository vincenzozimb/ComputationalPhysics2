#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "util.h"
#include "print_routines.h"

/* functions */
double potential(double r, double rho, double R){

    double pot;
    if(r <= R){
        pot = r * r / 3.0 - R * R;
    }else{
        pot = -2.0/3.0 * R * R * R / r; 
    }

    return 2.0 * M_PI * rho * pot;

}

double pot_eff(double r, double rho, double R, int l){
    return potential(r,rho,R) + (double)l*(l+1)/(2.0*r*r);
}

void print_potential(double L, double h, double rho, double R, int l){
    
    int dim = (int)(L/h);
    double r;

    FILE *file;
    file = fopen("potential.csv","w");
    
    for(int i=0; i<dim; i++){

        r = h * (i+1) + 0.5;
        fprint_double(file,r);
        fprint_double_newline(file,potential(r,rho,R));
    
    }
    
    fclose(file);
}

double F_pot(double r, void *p){

    Param *par = (Param *)p;

    double E = par->E;
    double rho = par->rho;
    double R = par->R;
    int l = par->l;

    return 2.0 * ( pot_eff(r,rho,R,l) - E );

}

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

double shooting(double E, void *p){

    /* extract parameters */
    Param *par = (Param*)p;
    par->E = E;
    double L = par->L;
    double h = par->h;
    int l = par->l;


    /* define variables */
    int dim = (int)(L/h);
    double r[dim], psi[dim];

    /* initial conditions */
    r[0] = h;
    r[1] = r[0] + h;

    psi[0] = pow(r[0],l+1);
    psi[1] = pow(r[1],l+1);

    /* evolve with numerov */
    initialize_position(r,h,dim);
    solve_numerov(r,psi,dim,h,F_pot,p);
    normalize(psi,h,dim);

    // FILE *file;
    // file = fopen("debug.csv","w");

    // fprint_two_vec(file,r,psi,dim);

    // fclose(file);

    /* return final value */
    return psi[dim-1];

}

double run_for_delta(double E, void *p){

    /* extract parameters */
    Param *par = (Param *)p;

    par->E = E;
    double L = par->L;
    double h = par->h;
    double rho = par->rho;
    double R = par->R;
    int l = par->l;


    /* calculate inversion point */
    // double x0 = zero_newton(F_pot,0.5,p);
    double x0 = sqrt(3) * sqrt(R*R + E/(2*M_PI*rho) );
    printf("x0 = %lf\n",x0);

    /* define variables */    
    int dimL = 1 + (int)floor(x0 / h);
    int dimR = 1 + (int)floor((L - x0) / h);

    double xL[dimL], xR[dimR], psiL[dimL], psiR[dimR];

    /* boundary conditions */
    xL[0] = h;
    xL[1] = xL[0] + h;

    xR[0] = L;
    xR[1] = xR[0] - h;

    psiL[0] = pow(xL[0],l+1);
    psiL[1] = pow(xL[1],l+1);

    psiR[0] = 1e230 * exp(2.0 * E * xR[0]);     // E is negative for bound states
    psiR[1] = 1e230 * exp(2.0 * E * xR[1]);
    printf("value = %lf\n",psiR[0]);

    /* evolution */
    initialize_position(xL,h,dimL);
    solve_numerov(xL,psiL,dimL,h,F_pot,p);

    initialize_position(xR,-h,dimR);
    solve_numerov(xR,psiR,dimR,-h,F_pot,p);

    /* normalize */
    normalize(psiL,h,dimL);
    normalize(psiR,h,dimR);

    debug_print(xR,psiR,dimR);
    debug_print(xL,psiL,dimL);

    /* adjust continuity */
    double K = psiL[dimL-1] / psiR[dimR-1];
    for(int i=0; i<dimR; i++){
        psiR[i] *= K;
    }

    /* return delta */
    assert( fabs(psiL[dimL-1] - psiR[dimR-1]) < EPS );
    double ris = psiL[dimL-2] + psiR[dimR-2];
    ris -= (2.0 + h * h * F_pot(x0,p)) * psiL[dimL-1];
    ris /= h;

    return ris;



}

void debug_print(double x[], double v[], int dim){
    FILE *file;
    file = fopen("debug.csv","w");
    fprint_two_vec(file,x,v,dim);
    fclose(file);
}

