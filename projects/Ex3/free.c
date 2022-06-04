/* ========================= ACTUALLY THIS IS THE NEW dft20.c ========================= */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <stdbool.h>

#include "lapack_wrappers.h"
#include "print_routines.h"

/* ======================= TYPEDEF ======================= */
typedef struct ParamPot{
    double R,rho;
}ParamPot;

typedef struct ParamEc{
    double p, A, alpha1, beta1, beta2, beta3, beta4;
}ParamEc;

/* ======================= GLOBAL VARIABLES ======================= */
#define N 20 // number of electrons
#define dim 500 // number of discretization points in the finite difference method
#define EPS 1e-4 // precision for the convergence of the autoconsistent cycle
#define BETA 1.0 // mixing procedure parameter

double rs, R, rho;
ParamPot par;

double L, h;
double r[dim];


/* ======================= FUNCTION HEADERS ======================= */
void fill_zero(double a[], int len);
void copy_vec(double copy[], double paste[], int len);
void fill_position();
double pot_ext(double r, void *p);
void add_pot_corr(double v[]);
void fill_pot_unch(double v[], bool free);
void add_pot_centr(double v[], int l);
void solve_radialSE_diagonalize(int Nb, double v[], double E[], double psi[][Nb]);
void normalize(int Nb, double psi[][Nb]);
void add_density(int Nb, double n[], double psi[][Nb], int l);
void solve_second_closed_shell(double E[4], double n[], double v[], bool free); // this is the only function specialized to the case N = 20
void density_integral(double n[]);
void add_pot_exc(double v[], double n[]);
void add_pot_coulomb(double v[], double n[]);
double L_one_distance(double a[], double b[], int len);
double spillout(double n[]);
void mixing_n(double n_old[], double n[]);

void print_func(double a[], double b[], int len, char name[25]);

/* ======================= MAIN ======================= */
int main(){

    /* system parameters */
    rs = 3.93; // 3.93 for Na and 4.86 for K
    R = rs * pow((double)N,1.0/3.0); // radius of the cluster (of its harmonic part)
    rho = 3.0 / (4.0 * M_PI * rs * rs * rs); // density of the jellium

    par.R = R;
    par.rho = rho;


    /* code parameters */
    L = 2.5 * R; // infrared cutoff. Space interval goes from 0 from L
    h = L / dim; // ultraviolet cutoff
   
    printf("\n=============== SYSTEM PARAMETERS ===============\n");
    printf("Number of electrons: N = %d\n",N);
    printf("The characteristic radius of the cluster is R = %lf\n",R);
    printf("The infrared cutoff is L = %lf\n",L);
    printf("The ultraviolet cutoff is h = %lf\n\n",h);

    bool free;


    /* variables */
    double v[dim], n_free[dim];
    double E[4];

    fill_position();


    /* Solving the SE for N = 20 non interacting electrons, in total we will have 4 orbitals (0s 0p 0d 1s) */
    fill_pot_unch(v,free=true);
    solve_second_closed_shell(E,n_free,v,free=true);
    print_func(r,n_free,dim,"density_free.csv");
    
    printf("=============== FREE ELECTRONS ===============\n");
    printf("The free energies are E_nl:\n");
    printf("E_%d%d = %lf\n",0,0,E[0]);
    printf("E_%d%d = %lf\n",1,0,E[1]);
    printf("E_%d%d = %lf\n",0,1,E[2]);
    printf("E_%d%d = %lf\n\n",0,2,E[3]);

    density_integral(n_free);


    /* Interacting electrons */
    double n_old[dim], n[dim], v_step[dim];
    double check;
    copy_vec(n_free,n_old,dim);
    fill_pot_unch(v,free=false);
    int cnt = 0;
    do{
        copy_vec(v,v_step,dim);
        add_pot_coulomb(v_step,n_old);
        add_pot_exc(v_step,n_old);
        copy_vec(n_old,n,dim);
        solve_second_closed_shell(E,n,v_step,free=false);
        check = L_one_distance(n,n_old,dim);
        copy_vec(n,n_old,dim);
        //mixing_n(n_old,n);
        cnt++;
        if(cnt == 5){
            // print_func(r,n,dim,"density.csv");
            // density_integral(n);
            break;
        }
        printf("\rcnt = %d, check = %lf",cnt,check);
        fflush(stdout);
    }while(check > EPS); 

    // print density and final energies
    printf("=============== INTERACTING ELECTRONS ===============\n");
    printf("The convercence was reached in %d steps\n",cnt);
    printf("The interacting eigenvalues are E_nl:\n");
    printf("E_%d%d = %lf\n",0,0,E[0]);
    printf("E_%d%d = %lf\n",1,0,E[1]);
    printf("E_%d%d = %lf\n",0,1,E[2]);
    printf("E_%d%d = %lf\n\n",0,2,E[3]);
    
    // print_func(r,n,dim,"density.csv");
    // density_integral(n);


    /* Calculate the spillout */
    double deltaN = spillout(n);
    printf("Spillout for N = 20: deltaN = %lf\n",deltaN);
    printf("\n");


}


/* ======================= FUNCTION BODIES ======================= */
// array functions
void fill_zero(double a[], int len){
    // initialize all the elements of the vector to zero
    for(int i=0; i<len; i++){
        a[i] = 0.0;
    }
}

void copy_vec(double copy[], double paste[], int len){
    // copy the vector copy[len] into the vector paste[len]
    for(int i=0; i<len; i++){
        paste[i] = copy[i];
    }
}

void fill_position(){
    // fill the position vector
    for(int i=0; i<dim; i++){
        r[i] = (i+1) * h;
    }
}

double pot_ext(double r, void *p){
    // analytic functional expression for the external potential
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

void add_pot_corr(double v[]){
    // add the correlation potential
    ParamEc p = {1.0, 0.031091, 0.21370, 7.5957, 3.5876, 1.6382, 0.49294};
    double K;
    K = 2.0 * p.A * (p.beta1 * pow(rs,0.5) + p.beta2 * rs + p.beta3 * pow(rs,1.5) + p.beta4 * pow(rs,p.p+1.0));
    K = 1.0 + 1.0/K; 
    K = log(K);
    K *= -2.0 * p.A * (1.0 + p.alpha1 * rs);

    for(int i=0; i<dim; i++){
        v[i] += K;
    }
}

void fill_pot_unch(double v[], bool free){
    // fill the vector v[dim] with the unchanged part of the potential
    for(int i=0; i<dim; i++){
        v[i] = pot_ext(r[i],&par);
    }
    if(!free){
        add_pot_corr(v);
    }
}

void add_pot_centr(double v[], int l){
    // add the centrifugal barrier to the potential vector
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

void normalize(int Nb, double psi[][Nb]){
    // Normalize the Nb orbitals contained in psi[dim][N]
    
    /* calculate norm */
    double norm[Nb];
    for(int i=0; i<Nb; i++){
        norm[i] = 0.0;
    }
    for(int i=0; i<Nb;i++){
        for(int j=0; j<dim; j++){
            norm[i] += (psi[j][i] * psi[j][i]) * pow(h*(j+1),2.0);
        }
        norm[i] *= h;
        norm[i] = sqrt(norm[i]);
    }
    /* normalize */
    for(int i=0; i<Nb; i++){
        for(int j=0; j<dim; j++){
            psi[j][i] /= norm[i];
        }
    }
}

void add_density(int Nb, double n[], double psi[][Nb], int l){
    // use the Nb orbitals in psi[dim][Nb] to update the density
    for(int i=0; i<dim; i++){
        for(int j=0; j<Nb; j++){
            n[i] += 2.0 * psi[i][j]*psi[i][j] * ((double)(2*l+1) / (4.0*M_PI));
        }
    }
}

void solve_second_closed_shell(double E[4], double n[], double v[], bool free){
    // diagonalize the SE for the (0s 0p 0d 1s) orbitals. Calculate the energies and fill the density array.
    double pot[dim];
    int l;
    
    int Nb;
    double *en, **psi;

    fill_zero(n,dim);

    // 0s and 1s
    copy_vec(v,pot,dim);
    l = 0; Nb = 2;

    en = malloc(Nb*sizeof(double));
    psi = malloc(dim*sizeof(double *));
    for(int i=0; i<dim; i++){
        psi[i] = malloc(Nb*sizeof(psi[i]));
    }

    add_pot_centr(pot,l);
    solve_radialSE_diagonalize(Nb,pot,en,psi);
    normalize(Nb,psi);
    add_density(Nb,n,psi,l);
    E[0] = en[0];
    E[1] = en[1];

    // 0p
    copy_vec(v,pot,dim);
    l = 1; Nb = 1;

    en = malloc(Nb*sizeof(double));
    psi = malloc(dim*sizeof(double *));
    for(int i=0; i<dim; i++){
        psi[i] = malloc(Nb*sizeof(psi[i]));
    }

    add_pot_centr(pot,l);
    solve_radialSE_diagonalize(Nb,pot,en,psi);
    normalize(Nb,psi);
    add_density(Nb,n,psi,l);
    E[2] = en[0];

    // 0d
    copy_vec(v,pot,dim);
    l = 2; Nb = 1;

    en = malloc(Nb*sizeof(double));
    psi = malloc(dim*sizeof(double *));
    for(int i=0; i<dim; i++){
        psi[i] = malloc(Nb*sizeof(psi[i]));
    }

    add_pot_centr(pot,l);
    solve_radialSE_diagonalize(Nb,pot,en,psi);
    normalize(Nb,psi);
    add_density(Nb,n,psi,l);
    E[3] = en[0];
   
}

void density_integral(double n[]){
    // Calculate the integral of the density and print the result. It should be equal to the number of electrons
    double tot = 0.0;
    for(int i=0; i<dim; i++){
        tot += n[i] * r[i] * r[i];
    }
    tot *= 4.0 * M_PI * h;
    printf("The integral of the density is: N = %lf\n\n",tot);
}

void add_pot_exc(double v[], double n[]){
    // add the exchange potential to the array v[dim]
    for(int i=0; i<dim; i++){
        v[i] += -3.0/4.0 * pow(3.0/M_PI,1.0/3.0) * pow(n[i],1.0/3.0);
    }
}

void add_pot_coulomb(double v[], double n[]){
    // add the Coulomb potential to the array v[dim]
    double K[dim];
    
    fill_zero(n,dim);
    n[100] = 1/h;

    print_func(r,n,dim,"density.csv");
    density_integral(n);

    fill_zero(K, dim);
    for(int i=0; i<dim; i++){
        for(int j=0; j<dim; j++){
            if(r[j]<r[i]){
                K[i] += r[j] * r[j] *n[j] / r[i];
            }else{
                K[i] += n[j] * r[j];
            }
        }
    }
    for(int i=0; i<dim; i++){
        K[i] *= 4.0 * M_PI * h;
    }
    print_func(r,K,dim,"data.csv");
    for(int i=0; i<dim; i++){
        v[i] += K[i];  
    }
}

double L_one_distance(double a[], double b[], int len){
    // Calculate the (discretized version of the) L1 distance between a[len] and b[len] 
    double dist = 0.0;
    for(int i=0; i<len; i++){
        dist += fabs(a[i]-b[i]);
    }
    dist *= h;
    dist = sqrt(dist);
    return dist;
}

double spillout(double n[]){
    // Calculate the electronic spillout associated to the given density n[dim]
    double ris = 0.0;
    for(int i=0; i<dim; i++){
        if(r[i]>R){
            ris += r[i] * r[i] * n[i];
        }
    }
    ris *= h * 4.0 * M_PI;
    return ris;
}

void mixing_n(double n_old[], double n[]){
    // Perform the mixing procedure
    for(int i=0; i<dim; i++){
        n_old[i] = BETA * n[i] + (1.0-BETA)*n_old[i];
    }
}

// print functions
void print_func(double a[], double b[], int len, char name[25]){
    // no input control
    // Print on file the two array
    FILE *file;
    file = fopen(name,"w");
    fprint_two_vec(file,a,b,len);
    fclose(file);
}
