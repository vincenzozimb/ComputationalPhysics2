#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "util.h"
#include "print_routines.h"
#include "lapack_wrappers.h"

/* global variables */
int N;
double rs, R, rho;

/* functions */
void solve_radialSE_diagonalize(int N, double r[], double v[], double E[], double psi[][N], double h, int dim){

    /* diagonal */
    double d[dim];
    for(int i=0; i<dim; i++){
        d[i] = 1/(h*h) + v[i];
    }

    /* subdiagonal */
    double sd[dim-1];
    for(int i=0; i<dim-1; i++){
        sd[i] = - 1.0 / (2.0 * h * h);
    }

    /* diagonalization */
    double eigval[dim];
    double eigvec[dim][dim];

    int info = diagonalize_tridiag_double(dim,d,sd,eigvec,eigval);
    assert(info == 0);

    /* save results */
    for(int i=0; i<N; i++){
        E[i] = eigval[i];
        for(int j=0; j<dim; j++){
            psi[j][i] = -eigvec[i][j] / r[j];   // in this way it return the R(r)
        }
    }

}

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

void normalize(int N, double psi[][N], double h, int dim){
    /* calculate norm */
    double norm[N];
    for(int i=0; i<N; i++){
        norm[i] = 0.0;
    }
    for(int i=0; i<N;i++){
        for(int j=0; j<dim; j++){
            norm[i] += (psi[j][i] * psi[j][i]) * pow(h*(j+1),2.0);
        }
        norm[i] *= h;
        norm[i] = sqrt(norm[i]);
    }
    /* normalize */
    for(int i=0; i<N; i++){
        for(int j=0; j<dim; j++){
            psi[j][i] /= norm[i];
        }
    }
}

void normalize_single(double psi[], double dx, int dim, double multiplier){
    /* calculate norm */
    double norm = 0.0;
    
    for(int i=0; i<dim; i++){
        norm += psi[i] * psi[i];
    }
    norm *= dx;
    norm = sqrt(norm);

    norm /= multiplier;

    /* normalize */
    for(int i=0; i<dim; i++){
        psi[i] /= norm;
    }
}

void print_func(double r[], double v[], int dim, char name[25]){

    // no controllo input

    FILE *file;
    file = fopen(name,"w");
    fprint_two_vec(file,r,v,dim);
    fclose(file);

}

void print_wf(int N, double r[], double psi[][N], int dim, double h, char name[6]){
    
    // no controllo input. Assume for name the spectroscopic notation

    FILE *file;
    file = fopen(name,"w");
    
    for(int i=0; i<dim; i++){
        for(int j=0; j<N+1; j++){
            if(j==0){
                fprint_double(file,r[i]);
            }else{
                fprint_double(file,psi[i][j-1]);
            }
        }
        fprintf(file,"\n");
    }

    fclose(file);

}

void fill_zero(double v[], int dim){
    for(int i=0; i<dim; i++){
        v[i] = 0.0;
    }
}

void fill_position(double r[], double h, int dim){
    
    for(int i=0; i<dim; i++){
        r[i] = h * (i+1);
    }
    
}

void fill_potential(double r[], double v[], int l, int dim, void *p){
        
    for(int i=0; i<dim; i++){
        v[i] = potential(r[i],p) + (double)l*(l+1)/(2.0*r[i]*r[i]);
    }
}

void add_density(int N, double r[], double n[], double psi[][N], int dim, int l){
    for(int i=0; i<dim; i++){
        for(int j=0; j<N; j++){
            n[i] += 2.0 * psi[i][j]*psi[i][j] * ((double)(2*l+1) / (4.0*M_PI));
        }
    }
}

void add_energy(int *cnt, double E[], double eps[], int Nb){

    for(int i=(*cnt); i<(*cnt)+Nb; i++){
        E[i] = eps[i-(*cnt)]; 
    }

}

void exchange_pot(double ex[], double n[], int dim){
    for(int i=0; i<dim; i++){
        ex[i] = -3.0/4.0 * pow(3.0/M_PI,1.0/3.0) * pow(n[i],1.0/3.0);
    }
}

void add_correlation_pot(double ec[], int dim){
    
    ParamEc p = {1.0, 0.031091, 0.21370, 7.5957, 3.5876, 1.6382, 0.49294};
    double K;
    K = 2.0 * p.A * (p.beta1 * pow(rs,0.5) + p.beta2 * rs + p.beta3 * pow(rs,1.5) + p.beta4 * pow(rs,p.p+1.0));
    K = 1.0 + 1.0/K; 
    K = log(K);
    K *= -2.0 * p.A * (1.0 + p.alpha1 * rs);

    for(int i=0; i<dim; i++){
        ec[i] += K;
    }
}

void add_coulomb_pot(double vh[], double r[], double n[], double h, int dim){

    double K = 0.0;
    for(int i=0; i<dim; i++){
        K += r[i] * r[i] * n[i];
    }
    K *= h;

    for(int i=0; i<dim; i++){
        vh[i] += K * 4.0 * M_PI / r[i];  
    } 

}

void partial_pot(double pot[], double r[], int dim, int l, void *p){

    // external and correlation potential, they do not depend on the density    
    fill_zero(pot,dim);
    add_correlation_pot(pot,dim);
    
    for(int i=0; i<dim; i++){
        pot[i] += potential(r[i],p);
    }

}

void change_pot(double pot[], double n[], double r[], double h, int dim){

    // the coulomb and the exchange potential do depend on the density
    add_coulomb_pot(pot,r,n,h,dim);
    add_correlation_pot(pot,dim);

}

void copy_vec(double copy[], double paste[], int dim){

    for(int i=0; i<dim; i++){
        paste[i] = copy[i];
    }

}

double euclidean_distance(double v1[], double v2[], int dim){

    double ris = 0.0;
    for(int i=0; i<dim; i++){
        ris += pow(v1[i]-v2[i],2.0);
    }
    ris = sqrt(ris);

    return ris;

}