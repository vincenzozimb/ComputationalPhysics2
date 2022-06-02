#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "util.h"
#include "print_routines.h"

/* ====================================================================== */
int main(){

    /* system parameters */
    int N = 40; // number of electrons
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
    double r[dim], v[dim], n_free[dim];

    fill_zero(n_free,dim);
    fill_position(r,h,dim);


    /* Solving the SE for N = 40 non interacting electrons, in total we will have 6 orbitals (0s 0p 0d 1s 0f 1p) */
    double E[6], eps[2];
    cnt = 0;
        
    for(int l=0; l<4; l++){
        if(l == 0){
            Nb = 2;
            char name[] = "s.csv"; 
            double psi[dim][Nb];
            fill_potential(r,v,l,dim,&par);
            solve_radialSE_diagonalize(Nb,r,v,eps,psi,h,dim);
            normalize(Nb,psi,h,dim);
            add_density(Nb,r,n_free,psi,dim,l);
            print_wf(Nb,r,psi,dim,h,name);
            add_energy(&cnt,E,eps,Nb);
            cnt += Nb;

        }else if(l==1){
            Nb = 2;
            char name[] = "p.csv";
            double psi[dim][Nb];
            fill_potential(r,v,l,dim,&par);
            solve_radialSE_diagonalize(Nb,r,v,eps,psi,h,dim);
            normalize(Nb,psi,h,dim);
            add_density(Nb,r,n_free,psi,dim,l);
            print_wf(Nb,r,psi,dim,h,name);
            add_energy(&cnt,E,eps,Nb);
            cnt += Nb;

        }else if(l==2){
            Nb = 1;
            char name[] = "d.csv";
            double psi[dim][Nb];
            fill_potential(r,v,l,dim,&par);
            solve_radialSE_diagonalize(Nb,r,v,eps,psi,h,dim);
            normalize(Nb,psi,h,dim);
            add_density(Nb,r,n_free,psi,dim,l);
            print_wf(Nb,r,psi,dim,h,name);
            add_energy(&cnt,E,eps,Nb);
            cnt += Nb;

        }else{
            Nb = 1;
            char name[] = "f.csv";
            double psi[dim][Nb];
            fill_potential(r,v,l,dim,&par);
            solve_radialSE_diagonalize(Nb,r,v,eps,psi,h,dim);
            normalize(Nb,psi,h,dim);
            add_density(Nb,r,n_free,psi,dim,l);
            print_wf(Nb,r,psi,dim,h,name);
            add_energy(&cnt,E,eps,Nb);
            cnt += Nb;
        }
    }

    printf("The energies E_ln are:\n");
    printf("E_%d%d = %lf\n",0,0,E[0]);
    printf("E_%d%d = %lf\n",0,1,E[1]);
    printf("E_%d%d = %lf\n",1,0,E[2]);
    printf("E_%d%d = %lf\n",1,1,E[3]);
    printf("E_%d%d = %lf\n",2,0,E[4]);
    printf("E_%d%d = %lf\n",3,0,E[5]);
    printf("\n");                                        


    /* print density and check its normalization */
    char name[] = "density_free.csv";
    print_func(r,n_free,dim,name);

    double tot = 0.0;
    for(int i=0; i<dim; i++){
        tot += n_free[i] * r[i] * r[i];
    }
    tot *= 4.0 * M_PI * h;
    printf("The intregral of the density is: N = %lf\n",tot);

    printf("\n");


}