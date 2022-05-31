#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "util.h"
#include "print_routines.h"

/* ====================================================================== */
int main(){

    /* system parameters */
    int N = 8; // number of electrons
    double rs = 3.93; // 3.93 for Na and 4.86 for K

    double R = rs * pow((double)N,1.0/3.0); // radius of the cluster (of its harmonic part)
    double rho = 3.0 / (4.0 * M_PI * rs * rs * rs); // density of the jellium

    ParamPot par = {R,rho};


    /* code parameters */
    double L = 2.0 * R; // infrared cutoff. Space interval goes from 0 from L
    int dim = 900; // number of discretizations in the finite difference method
    double h = L / dim; // ultraviolet cutoff

    printf("\nThe characteristic radius of the cluster is R = %lf\n",R);
    printf("The infrared cutoff is L = %lf\n",L);
    printf("The ultraviolet cutoff is h = %lf\n\n",h);

    int Nb, cnt;
    double r[dim], v[dim], n[dim];

    for(int i=0; i<dim; i++){
        n[i] = 0.0;
    }

    fill_position(r,h,dim);

    
    /* Solving the SE for N = 8 non interacting electrons, in total we will have 2 orbitals (0s 0p) */
    double E[2], eps[1];
    cnt = 0;
    
    for(int l=0; l<2; l++){
        if(l == 0){
            Nb = 1;
            char name[] = "s.csv"; 
            double psi[dim][Nb];
            fill_potential(r,v,l,dim,&par);
            solve_radialSE_diagonalize(Nb,r,v,eps,psi,h,dim);
            normalize(Nb,psi,h,dim);
            add_density(Nb,r,n,psi,dim,l);
            print_wf(Nb,r,psi,dim,h,name);
            add_energy(&cnt,E,eps,Nb);
            cnt += Nb;

        }else{
            Nb = 1;
            char name[] = "p.csv";
            double psi[dim][Nb];
            fill_potential(r,v,l,dim,&par);
            solve_radialSE_diagonalize(Nb,r,v,eps,psi,h,dim);
            normalize(Nb,psi,h,dim);
            add_density(Nb,r,n,psi,dim,l);
            print_wf(Nb,r,psi,dim,h,name);
            add_energy(&cnt,E,eps,Nb);
            cnt += Nb;

        }
    }

    printf("The energies E_ln are:\n");
    printf("E_%d%d = %lf\n",0,0,E[0]);
    printf("E_%d%d = %lf\n",1,0,E[1]);
    printf("\n"); 


    /* print density and check its normalization */
    char name[] = "density.csv";
    print_func(r,n,dim,name);

    double tot = 0.0;
    for(int i=0; i<dim; i++){
        tot += n[i] * r[i] * r[i];
    }
    tot *= 4.0 * M_PI * h;
    printf("The intregral of the density is: N = %lf\n",tot);

    printf("\n");


}