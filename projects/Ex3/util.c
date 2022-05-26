#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "util.h"
#include "print_routines.h"
#include "lapack_wrappers.h"


/* functions */
void solve_radialSE_diagonalize(int N, int l, double E[], double psi[][N], double L, int dim, double pot(double, void*), void *p, int bol){

    double h = L / dim; // ultraviolet cutoff

    double r[dim], v[dim], d[dim], sd[dim-1];

    /* position and potential */
    
    for(int i=0; i<dim; i++){
        r[i] = h * (i+1);
        v[i] = pot(r[i],p) + (double)l*(l+1)/(2.0*r[i]*r[i]);
    }

    /* print potential */
    if(bol == 1){
        FILE *file;
        file = fopen("potential.csv","w");
        fprint_two_vec(file,r,v,dim);
        fclose(file);
    }

    /* diagonal */
    for(int i=0; i<dim; i++){
        d[i] = 1/(h*h) + v[i];
    }

    /* subdiagonal */
    for(int i=0; i<dim-1; i++){
        sd[i] = - 1.0 / (2.0 * h * h);
    }

    /* diagonalization */
    double eigval[dim];
    double eigvec[dim][dim];

    int info = diagonalize_tridiag_double(dim,d,sd,eigvec,eigval);
    assert(info == 0);

    /* print result */
    int cnt = 0;
    for(int i=0; i<dim; i++){
        if(cnt < N && eigval[i] < 0.0){
            printf("E_nl \t E_%d%d = %.3lf\n",cnt,l,eigval[i]);
            E[cnt] = eigval[i];
            cnt++;
        }
    }

    for(int j=0; j<N; j++){
        for(int i=0; i<dim; i++){
            psi[i][j] = -eigvec[j][i];
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

double harmonic(double r, void *p){

    return 0.5 * r * r - 5.0;

}

void normalize(int N, double psi[][N], double h, int dim){
    
    /* calculate norm */
    double norm[N];
    for(int i=0; i<N; i++){
        norm[i] = 0.0;
    }
    
    for(int i=0; i<N;i++){
        for(int j=0; j<dim; j++){
            norm[i] += psi[j][i] * psi[j][i];
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

void normalize_single(double psi[], double dx, int dim){
    /* calculate norm */
    double norm = 0.0;
    
    for(int i=0; i<dim; i++){
        norm += psi[i] * psi[i];
    }
    norm *= dx;
    norm = sqrt(norm);

    norm /= 10.0;

    /* normalize */
    for(int i=0; i<dim; i++){
        psi[i] /= norm;
    }
}

void print_wf(int N, double psi[][N], int dim, double h){
    
    FILE *file;
    file = fopen("wf.csv","w");
    
    for(int i=0; i<dim; i++){
        for(int j=0; j<N+1; j++){
            if(j==0){
                fprint_double(file,h*(i+1));
            }else{
                fprint_double(file,psi[i][j]);
            }
        }
        fprintf(file,"\n");
    }

    fclose(file);

}

