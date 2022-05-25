#include <stdio.h>
#include <math.h>
#include <assert.h>

#include "lapack_wrappers.h"
#include "print_routines.h"

/* =============== */
// double potential(double r, int l);
double potential(double r, double rho, double R, int l);
double K(double x, double rho, double R, double h, int l);
void print_potential(double x[], double v[], int dim);
void print_wf(double x[], double v[], int dim);
void normalize(double psi[], double dx, int dim);

/* =============== */
int main(){
 
    /* parameters */
    int N = 20; // number of electrons
    double rs = 3.93; // 3.93 for Na and 4.86 for K

    double R = rs * pow((double)N,1.0/3.0); // radius of the cluster (of its harmonic part)
    double rho = 3.0 / (4.0 * M_PI * rs * rs * rs); // density of the jellium

    double L = 3.0 * R;
    int dim = 950;
    double h = L / dim;

    printf("\nThe ultraviolet cutoff is h=%lf\n",h);
    printf("R=%lf\n",R);
    printf("L=%lf\n\n",L);

    double r0 = 0.0; // r0 != 0 only for plotting purposes

    double r[dim], v[dim], d[dim], sd[dim-1];

    int l = 0;

    /* position and potential */
    for(int i=0; i<dim; i++){
        r[i] = r0 + h * (i+1);
        v[i] = potential(r[i],rho,R,l);
    }

    print_potential(r,v,dim);

    /* diagonal */
    for(int i=0; i<dim; i++){
        d[i] = K(r[i],rho,R,h,l);
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
        if(cnt < 6){
            printf("l=%d \t E_%d = %.3lf\n",l,cnt,eigval[i]);
            cnt++;
        }
    }

    double psi[dim];
    for(int i=0; i<dim; i++){
        psi[i] = -eigvec[0][i];
    }
    normalize(psi,h,dim);
    print_wf(r,psi,dim);

}


/* =============== */
// double potential(double r, int l){
//     return r * r / 2.0 + (double)l*(l+1)/(2*r*r);
// }

double potential(double r, double rho, double R, int l){

    double pot;
    if(r <= R){
        pot = r * r / 3.0 - R * R;
    }else{
        pot = -2.0/3.0 * R * R * R / r; 
    }

    pot *= 2.0 * M_PI * rho;
    pot += (double)l*(l+1) / (2.0*r*r);
    return pot;

}

double K(double x, double rho, double R, double h, int l){
    return potential(x,rho,R,l) + 1 / (h * h);
}

void print_potential(double x[], double v[], int dim){

    FILE *file;
    file = fopen("potential.csv","w");
    fprint_two_vec(file,x,v,dim);
    fclose(file);

}

void print_wf(double x[], double v[], int dim){

    FILE *file;
    file = fopen("wf.csv","w");
    fprint_two_vec(file,x,v,dim);
    fclose(file);

}

void normalize(double psi[], double dx, int dim){
    /* calculate norm */
    double norm = 0.0;
    
    for(int i=0; i<dim; i++){
        norm += psi[i] * psi[i];
    }
    norm *= dx;
    norm = sqrt(norm);

    norm *= 1.0 / 10.0;

    /* normalize */
    for(int i=0; i<dim; i++){
        psi[i] /= norm;
    }
}