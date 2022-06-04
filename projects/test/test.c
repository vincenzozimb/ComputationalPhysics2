#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <stdbool.h>

#include "lapack_wrappers.h"
#include "print_routines.h"

/* ======================= TYPEDEF ======================= */

/* ======================= GLOBAL VARIABLES ======================= */
#define dim 500 // discretization

double L,h;
double r[dim];

/* ======================= FUNCTION HEADERS ======================= */
double pot(double r);
void fill_position();
void fill_potential(double v[]);
void add_centr(double v[], int l);
void solve_radialSE_diagonalize(int Nb, double v[], double E[], double psi[][Nb]);
void print_func(double a[], double b[], int len, char name[25]);
void print_wf(int Nb, double psi[][Nb], char name[25]);

/* ======================= MAIN ======================= */
int main(){

    /* parameters */
    L = 8.0; // infrared cutoff
    h = (L/dim); // ultraviolet cutoff

    printf("\nThe infrared cutoff is L = %lf\n",L);
    printf("The ultraviolet cutoff is h = %lf\n\n",h);


    /* variables */
    int Nb = 4;
    double E[Nb];
    double v[dim], psi[dim][Nb];

    fill_position();
    fill_potential(v);
    
    int l = 0;
    add_centr(v,l);

    solve_radialSE_diagonalize(Nb,v,E,psi);

    printf("\nThe energies E_ln are:\n");
    printf("E_%d%d = %lf\n",0,0,E[0]);
    printf("E_%d%d = %lf\n",0,1,E[1]);
    printf("E_%d%d = %lf\n",0,2,E[2]);
    printf("E_%d%d = %lf\n",0,3,E[3]);
    printf("\n");  

    print_wf(Nb,psi,"data.csv");    

}


/* ======================= FUNCTION BODIES ======================= */
double pot(double r){
    return 0.5 * r * r;
}

void fill_position(){
    for(int i=0; i<dim; i++){
        r[i] = h * (i+1);
    }
}

void fill_potential(double v[]){
    for(int i=0; i<dim; i++){
        v[i] = pot(r[i]);
    }
}

void add_centr(double v[], int l){
    for(int i=0; i<dim; i++){
        v[i] += (double)l*(l+1)/(2.0*r[i]*r[i]);
    }
}

void solve_radialSE_diagonalize(int Nb, double v[], double E[], double psi[][Nb]){
    // Diagonalize the SE with the finite difference method (i.e. using as a basis the position eigenstates),
    // finding the first Nb bound states of the potential in v[dim], and save the results in E[dim] and psi[dim][Nb]

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
    for(int i=0; i<Nb; i++){
        E[i] = eigval[i];
        for(int j=0; j<dim; j++){
            psi[j][i] = -eigvec[i][j] / r[j];   // in this way it return the R(r)
        }
    }

}

void print_func(double a[], double b[], int len, char name[25]){
    // no input control
    // Print on file the two array
    FILE *file;
    file = fopen(name,"w");
    fprint_two_vec(file,a,b,len);
    fclose(file);
}

void print_wf(int Nb, double psi[][Nb], char name[25]){
    FILE *file;
    file = fopen(name,"w");
    for(int i=0; i<dim; i++){
        for(int j=0; j<Nb+1; j++){
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