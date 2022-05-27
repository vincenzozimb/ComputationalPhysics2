#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "util.h"
#include "print_routines.h"

/* ====================================================================== */
int main(){

    /* system parameters */
    int N = 20; // number of electrons
    double rs = 3.93; // 3.93 for Na and 4.86 for K

    double R = rs * pow((double)N,1.0/3.0); // radius of the cluster (of its harmonic part)
    double rho = 3.0 / (4.0 * M_PI * rs * rs * rs); // density of the jellium

    ParamPot par = {R,rho};


    /* code parameters */
    double L = 3.0 * R; // infrared cutoff. Space interval goes from 0 from L
    int dim = 900; // number of discretizations in the finite difference method
    double h = L / dim; // ultraviolet cutoff

    printf("\nThe characteristic radius of the cluster is R = %lf\n",R);
    printf("The infrared cutoff is L = %lf\n",L);
    printf("The ultraviolet cutoff is h = %lf\n\n",h);


    /* Solving the SE for the non interacting electrons, in total we will have 4 orbitals (0s 0p 0d 1s) */
    double r[dim], E[4], psi[dim][4];
    fill_position(r,h,dim);

    /* 0s and 1s orbitals */
    int ls = 0;
    int Nbs = 2;

    double vs[dim], Es[Nbs], psis[dim][Nbs];

    fill_potential(r,vs,ls,dim,&par);
    solve_radialSE_diagonalize(Nbs,r,vs,Es,psis,h,dim);    
    normalize(Nbs,psis,h,dim);
    fprint_vec(stdout,Es,Nbs);

    for(int i=0; i<Nbs; i++){
        E[i] = Es[i];
    }
    for(int i=0; i<dim; i++){
        for(int j=0; j<Nbs; j++){
            psi[i][j] = psis[i][j];
        }
    }

    /* 0p orbital */
    int lp = 1;
    int Nbp = 1;

    double vp[dim], Ep[Nbp], psip[dim][Nbp];

    fill_potential(r,vp,lp,dim,&par);
    solve_radialSE_diagonalize(Nbp,r,vp,Ep,psip,h,dim);    
    normalize(Nbp,psip,h,dim);
    fprint_vec(stdout,Ep,Nbp);

    for(int i=Nbs; i<Nbs+Nbp; i++){
        E[i] = Ep[i];
    }
    for(int i=0; i<dim; i++){
        for(int j=Nbs; j<Nbs+Nbp; j++){
            psi[i][j] = psip[i][j];
        }
    }

    /* 0d orbital */
    int ld = 2;
    int Nbd = 1;

    double vd[dim], Ed[Nbd], psid[dim][Nbd];

    fill_potential(r,vd,ld,dim,&par);
    solve_radialSE_diagonalize(Nbd,r,vd,Ed,psid,h,dim);    
    normalize(Nbd,psid,h,dim);
    fprint_vec(stdout,Ed,Nbd);
    
    for(int i=Nbs+Nbp; i<Nbs+Nbp+Nbd; i++){
        E[i] = Ed[i];
    }
    for(int i=0; i<dim; i++){
        for(int j=Nbs+Nbp; j<Nbs+Nbp+Nbd; j++){
            psi[i][j] = psid[i][j];
        }
    }

    print_wf(4,r,psi,dim,h);
    fprint_vec(stdout,E,4);


    /* density */
    // double n[dim];

    // density(Nb,r,n,psi,h,dim);
    // FILE *file;
    // file = fopen("density.csv","w");
    // fprint_two_vec(file,r,n,dim);
    // fclose(file);

    // double norm = 0.0;
    // for(int i=0; i<dim; i++){
    //     norm += n[i]*n[i];
    // }
    // norm *= h;
    // printf("Norm = %lf\n",norm);


    printf("\n");

}

/* ====================================================================== */
