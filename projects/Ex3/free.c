#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "print_routines.h"
#include "lapack_wrappers.h"

/*=============== STRUCT ===============*/
typedef struct Param_pot{
    double rho, R;
    int l;
}Param_pot;


/*=============== FUNCTION HEADERS ===============*/
double potential(double r, Param_pot *p);
double pot_eff(double r, Param_pot *p);
void fill_vec(double r[], double v[], double h, int dim, Param_pot *p);
double K(double r, double h, Param_pot *p);
void print_pot(double r[], double v[], int dim);

/*=============== MAIN ===============*/
int main(){

    /* system parameters */
    int N = 20; // number of electrons
    double rs = 3.93; // 3.93 for Na and 4.86 for K

    double R = rs * pow((double)N,1.0/3.0); // radius of the cluster (of its harmonic part)
    double rho = 3.0 / (4.0 * M_PI * rs * rs * rs); // density of the jellium

    /* code parameters */
    double L = 3.0 * R; // infrared cutoff. Space interval goes from 0 from L
    int dim = 500;

    double h = L / dim; // ultraviolet cutoff

    printf("R = %lf\n",R);
    printf("L = %lf\n",L);

    double r[dim], v[dim], diag[dim], subdiag[dim-1], psi[dim];

    int l = 0;
    Param_pot par = {rho,R,l};

    /* fill vectors */
    fill_vec(r,v,h,dim,&par);
    print_pot(r,v,dim);

    /* initialize matrix */
    for(int i=0; i<dim-1; i++){
        subdiag[i] = -1.0;
        diag[i] = K(r[i],h,&par); 
    }
    diag[dim-1] = K(r[dim-1],h,&par);

    /* diagonalization */
    double eigvec[dim][dim];
    double eigval[dim];

    int info = diagonalize_tridiag_double(dim,diag,subdiag,eigvec,eigval);
    assert(info==0);

    /* print result */
    int cnt = 0;
    // int nbs = 3;
    for(int i=0;i<dim;i++){
        if(eigval[i]<0.0){
            printf("E_%d = %lf\n",cnt,eigval[i]/(2.0*h*h));
            cnt++;
        }
    }





}


/*=============== FUNCTION BODIES ===============*/
// double potential(double r, Param_pot *p){
//     return 0.5 * r * r;
// }

double potential(double r, Param_pot *p){

    double rho = p->rho;
    double R = p->R;

    double pot;
    if(r <= R){
        pot = r * r / 3.0 - R * R;
    }else{
        pot = -2.0/3.0 * R * R * R / r; 
    }

    return 2.0 * M_PI * rho * pot;

}

void print_pot(double r[], double v[], int dim){
    FILE *file;
    file = fopen("potential.csv","w");
    fprint_two_vec(file,r,v,dim);
    fclose(file);
}
double pot_eff(double r, Param_pot *p){
    
    int l = p->l;

    return potential(r,p) + (double)l*(l+1)/(2.0*r*r);
}

void fill_vec(double r[], double v[], double h, int dim, Param_pot *p){
    
    for(int i=0; i<dim; i++){
        r[i] = h * (i + 1);
        v[i] = pot_eff(r[i],p);
    }

}

double K(double r, double h, Param_pot *p){

    return 2.0 + 2.0 * h * h * pot_eff(r,p);

}