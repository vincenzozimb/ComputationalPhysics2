#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "print_routines.h"

/* ======================= TYPEDEF ======================= */
typedef struct ParamPot{
    double R,rho;
}ParamPot;


/* ======================= GLOBAL VARIABLES ======================= */
int N;
double rs, R, rho;
ParamPot par;


/* ======================= FUNCTION HEADERS ======================= */
void fill_position(double r[], double h, int dim);
double potential(double r);
void fill_potential(double v[], double r[], int l, int dim);
void fill_F(double F[], double v[], double E, int dim);
void initial_conditions(double r[], double psi[], int l);
void fill_psi(double psi[], double F[], double h, int dim);
void normalize(double psi[], double h, int dim);
double shooting(double E, double r[], double v[], double h, int l, int dim);


/* ======================= MAIN ======================= */
int main(){

    /* system parameters */
    N = 20; // number of electrons
    rs = 3.93; // 3.93 for Na and 4.86 for K

    R = rs * pow((double)N,1.0/3.0); // radius of the cluster (of its harmonic part)
    rho = 3.0 / (4.0 * M_PI * rs * rs * rs); // density of the jellium

    par.R = R;
    par.rho = rho;


    /* code parameters */
    double L = 2.5 * R; // infrared cutoff. Space interval goes from 0 from L
    double h = 1e-3; // ultraviolet cutoff
    int dim = (int)(L/h);
   
    printf("\nThe characteristic radius of the cluster is R = %lf\n",R);
    printf("The infrared cutoff is L = %lf\n",L);
    printf("The ultraviolet cutoff is h = %lf\n\n",h);


    /* variables */
    double r[dim], v[dim];

    int l = 0;

    fill_position(r,h,dim);
    fill_potential(v,r,l,dim);
    
    /* solving */
    double E = -2.8;
    double F[dim], psi[dim];

    FILE *file;
    file = fopen("wf.csv","w");

    fill_F(F,v,E,dim);
    initial_conditions(r,psi,l);
    fill_psi(psi,F,h,dim);
    normalize(psi,h,dim);
    fprint_two_vec(file,r,psi,dim);

    fclose(file);

    double dE = 0.01;
    E = -2.0 * M_PI * rho * R * R + dE;
    double E_max = 0.0;

    file = fopen("data.csv","w");
    do{
        fprint_double(file,E);
        fprint_double_newline(file,shooting(E,r,v,h,l,dim));
        printf("E = %lf\n",E);
        E += dE;
    }while(E<E_max);
    fclose(file);












    // FILE *file;
    // file = fopen("data.csv","w");
    // fprint_two_vec(file,r,v,dim);
    // fclose(file);

}


/* ======================= FUNCTION BODIES ======================= */
void fill_position(double r[], double h, int dim){
    for(int i=0; i<dim; i++){
        r[i] = (i+1) * h;
    }
}

double potential(double r){

    double pot;
    if(r <= R){
        pot = r * r / 3.0 - R * R;
    }else{
        pot = -2.0/3.0 * R * R * R / r; 
    }

    pot *= 2.0 * M_PI * rho;
    return pot;

}

void fill_potential(double v[], double r[], int l, int dim){
    for(int i=0; i<dim; i++){
        v[i] = potential(r[i]) + (double)l*(l+1) / (2.0*r[i]*r[i]);
    }
}

void fill_F(double F[], double v[], double E, int dim){
    for(int i=0; i<dim; i++){
        F[i] = 2.0 * (v[i] - E);
    }
}

void initial_conditions(double r[], double psi[], int l){
    psi[0] = pow(r[0],l+1);
    psi[1] = pow(r[1],l+1);
}

void fill_psi(double psi[], double F[], double h, int dim){
    for(int i=2; i<dim; i++){
        psi[i] = (2.0 + 5.0/6.0*h*h*F[i-1])*psi[i-1];
        psi[i] -= (1.0 -h*h/12.0*F[i-2])*psi[i-2];
        psi[i] /= 1.0-h*h/12.0*F[i];
    }
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

double shooting(double E, double r[], double v[], double h, int l, int dim){

    double F[dim], psi[dim];

    fill_F(F,v,E,dim);
    initial_conditions(r,psi,l);
    fill_psi(psi,F,h,dim);
    normalize(psi,h,dim);

    return psi[dim-1];

}

