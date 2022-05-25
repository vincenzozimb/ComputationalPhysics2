#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "util.h"
#include "print_routines.h"


int main(){

    /* system parameters */
    int N = 20; // number of electrons
    double rs = 3.93; // 3.93 for Na and 4.86 for K

    double R = rs * sqrt((double)N); // radius of the cluster (of its harmonic part)
    double rho = 3.0 / (4.0 * M_PI * rs * rs * rs); // density of the jellium


    /* code parameters */
    double L = 2.5 * R; // infrared cutoff. Space interval goes from 0 from L
    double h = 1e-3; // ultraviolet cutoff

    printf("R = %lf\n",R);
    printf("L = %lf\n",L);
    printf("depth for l=0 = %lf\n",2*M_PI*rho*R*R);


    /* number of bound state to be calculed */
    int l = 1;

    print_potential(L,h,rho,R,l);

    double E = -6.0;
    double E_max = -h; // it has to be less than the depth of the potential well
    double dE = 0.005;

    // Param par; //= {L,h,E,xi,R,rho,l};
    // par.L = L;
    // par.h = h;
    // par.E = E;
    // par.R = R;
    // par.rho = rho;
    // par.l = l;

    //double delta = run_for_delta(E,&par);

    // FILE *file;
    // file = fopen("debug.csv","w");

    // int dim = (int)((E_max-E)/dE);
    // for(int i=0; i<dim; i++){
    //     fprint_double(file,E);
    //     fprint_double_newline(file,run_for_delta(E,&par));
    //     E += dE;
    // }

    // fclose(file);







}