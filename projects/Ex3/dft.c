#include <stdio.h>
#include <stdlib.h>
#include <math.h>

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

    printf("\nThe ultraviolet cutoff is h = %lf\n",h);
    printf("R = %lf\n",R);
    printf("L = %lf\n\n",L);

    int l = 0;  // angular momentum
    int Nb = 3;  // number of bound state to be found
    
    double E[Nb], psi[dim][Nb];

    solve_radialSE_diagonalize(Nb,l,E,psi,L,dim,potential,&par,1);
        
    normalize(Nb,psi,h,dim);
    print_wf(Nb,psi,dim,h);

    printf("\n");

}

/* ====================================================================== */
