#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "util.h"
#include "print_routines.h"
#include "bisection.h"

int main(){

    /* parameters */
    double L = 6.0; // space interval goes from -L to L
    double h = 1e-3;

    /* number of bound state to be calculed */
    int n = 5;
    
    double E = 0.1;
    double dE = 0.05;
    double E_bound[n];

    /* initialization */
    Param_delta pd;
    pd.L = L;
    pd.h = h;
    pd.n = 0; // start with the ground state
    
    /* run and save */    
    double delta_min, delta, delta_new;
    delta_min = run_for_delta(E,&pd);

    do{

        delta = run_for_delta(E,&pd);
        delta_new = run_for_delta(E+dE,&pd);

        if(delta_min * delta_new < 0.0){
            E_bound[pd.n] = bisection(run_for_delta,E-dE,E+dE,&pd);
            pd.n++;
            delta_min = delta;
        }

        E += dE;

    }while(pd.n < n);

    printf("\nThe eigenvalues are:\n");
    for(int i=0; i<n; i++){
        printf("E_%d = %lf\n",i,E_bound[i]);
    }
    
    /* eigenfunctions */
    int dim = (2.0 * L / h);
    double x[dim], psi[dim];

    printf("\nWhich eigenfuncion do you want to plot? (choose a number between 0 and 4) \nn = ");
    int nb; // eigenfunction to be plotted
    scanf("%d",&nb);
    assert(nb < n);

    Param_F pf;

    initial_condition_ho(x,psi,L,h,nb);
    initialize_position(x,h,dim);

    pf.E = E_bound[nb];
    solve_numerov(x,psi,dim,h,F_ho,&pf);
    normalize(psi,h,dim);

    /* save data */
    FILE *file;
    file = fopen("psi1d.csv","w");
    fprint_two_vec(file,x,psi,dim);
    fclose(file);

}
