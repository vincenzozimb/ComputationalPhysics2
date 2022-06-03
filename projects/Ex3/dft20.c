#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "util.h"
#include "print_routines.h"

/* ====================================================================== */
int main(){

    /* system parameters */
    N = 20; // number of electrons
    rs = 3.93; // 3.93 for Na and 4.86 for K

    R = rs * pow((double)N,1.0/3.0); // radius of the cluster (of its harmonic part)
    rho = 3.0 / (4.0 * M_PI * rs * rs * rs); // density of the jellium

    ParamPot par = {R,rho};


    /* code parameters */
    double L = 2.0 * R; // infrared cutoff. Space interval goes from 0 from L
    int dim = 200; // number of discretizations in the finite difference method
    double h = L / dim; // ultraviolet cutoff

    printf("\nThe characteristic radius of the cluster is R = %lf\n",R);
    printf("The infrared cutoff is L = %lf\n",L);
    printf("The ultraviolet cutoff is h = %lf\n\n",h);

    int Nb, cnt;
    double r[dim], v[dim], n_free[dim];

    fill_zero(n_free,dim);
    fill_position(r,h,dim);

    /* Solving the SE for N = 20 non interacting electrons, in total we will have 4 orbitals (0s 0p 0d 1s) */
    double E[4], eps[2];
    cnt = 0;

    for(int l=0; l<3; l++){
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
            Nb = 1;
            char name[] = "p.csv";
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
            char name[] = "d.csv";
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
    printf("E_%d%d = %lf\n",2,0,E[3]);
    printf("\n");  
    
    /* print free density and check its normalization */
    char name[] = "density_free.csv";
    print_func(r,n_free,dim,name);

    double tot = 0.0;
    for(int i=0; i<dim; i++){
        tot += n_free[i] * r[i] * r[i];
    }
    tot *= 4.0 * M_PI * h;
    printf("The intregral of the density is: N = %lf\n",tot);

    printf("\n");


    /* Interacting electrons */
    double n[dim];
    double E_new[4];
    fill_zero(E_new,4);

    copy_vec(n_free,n,dim);

    double check, Vs[dim], Vp[dim], Vd[dim];
    partial_pot(Vs,r,dim,0,&par);
    partial_pot(Vp,r,dim,1,&par);
    partial_pot(Vd,r,dim,2,&par);

    double vs[dim], vp[dim], vd[dim];

    do{

        fill_zero(E,4);
        fill_zero(eps,2);
        cnt = 0;

        for(int l=0; l<3; l++){
            if(l == 0){
                Nb = 2;
                double psi[dim][Nb];
                copy_vec(Vs,vs,dim);
                change_pot(vs,n,r,h,dim);
                solve_radialSE_diagonalize(Nb,r,vs,eps,psi,h,dim);            
                normalize(Nb,psi,h,dim);
                add_density(Nb,r,n,psi,dim,l);
                add_energy(&cnt,E,eps,Nb);
                cnt += Nb;

            }else if(l==1){
                Nb = 1;
                double psi[dim][Nb];
                copy_vec(Vp,vp,dim);
                change_pot(vp,n,r,h,dim);
                solve_radialSE_diagonalize(Nb,r,vp,eps,psi,h,dim);
                normalize(Nb,psi,h,dim);
                add_density(Nb,r,n,psi,dim,l);
                add_energy(&cnt,E,eps,Nb);
                cnt += Nb;

            }else{
                Nb = 1;
                double psi[dim][Nb];
                copy_vec(Vd,vd,dim);
                change_pot(vd,n,r,h,dim);
                solve_radialSE_diagonalize(Nb,r,vd,eps,psi,h,dim);
                normalize(Nb,psi,h,dim);
                add_density(Nb,r,n,psi,dim,l);
                add_energy(&cnt,E,eps,Nb);
                cnt += Nb;

            }
        }

        check = euclidean_distance(E_new,E,4);
        copy_vec(E,E_new,4);

        printf("dist = %lf\n",check);

    }while(check > EPS);



}