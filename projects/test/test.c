#include <stdio.h>
#include <math.h>
#include <assert.h>

#include "lapack_wrappers.h"
#include "print_routines.h"

/* =============== */
typedef struct ParamPot{
    double R, rho;
}ParamPot;

/* =============== */
void solve_radialSE_diagonalize(int N, double r[], double v[], double E[], double psi[][N], double h, int dim);
double potential(double r, void *p);
void fill_position(double r[], double h, int dim);
void fill_potential(double r[], double v[], int l, int dim, void *p);
void print_potential(double x[], double v[], int dim);
void print_wf(int N, double r[], double psi[][N], int dim, double h);
void normalize(int N, double psi[][N], double h, int dim);

/* =============== */
int main(){
 
    /* parameters */
    int N = 20; // number of electrons
    double rs = 3.93; // 3.93 for Na and 4.86 for K

    double R = rs * pow((double)N,1.0/3.0); // radius of the cluster (of its harmonic part)
    double rho = 3.0 / (4.0 * M_PI * rs * rs * rs); // density of the jellium

    ParamPot par = {R,rho};

    double L = 3.0 * R; // infrared cutoff
    int dim = 950;
    double h = L / dim; // ultraviolet cutoff

    printf("\nThe ultraviolet cutoff is h=%lf\n",h);
    printf("R=%lf\n",R);
    printf("L=%lf\n\n",L);

    /* Solving the SE for the non interacting electrons, in total we will have 4 orbitals (0s 0p 0d 1s) */
    int l = 2;
    int Nb = 1;
    double r[dim], E[Nb], psi[dim][Nb],v[dim];

    fill_position(r,h,dim);
    fill_potential(r,v,l,dim,&par);

    solve_radialSE_diagonalize(Nb,r,v,E,psi,h,dim);    
    normalize(Nb,psi,h,dim);
    print_wf(Nb,r,psi,dim,h);

   
}

/* =============== */
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
            psi[j][i] = -eigvec[i][j];
        }
    }

}

double potential(double r, void *p){

    ParamPot *par = (ParamPot*)p;
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

void print_potential(double x[], double v[], int dim){

    FILE *file;
    file = fopen("potential.csv","w");
    fprint_two_vec(file,x,v,dim);
    fclose(file);

}

void print_wf(int N, double r[], double psi[][N], int dim, double h){
    
    FILE *file;
    file = fopen("wf.csv","w");
    
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
