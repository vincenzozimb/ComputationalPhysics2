#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "util.h"
#include "print_routines.h"
#include "bisection.h"

int main(){
    /* parameters */
    double L = 6.0;
    double h = 1e-3;
    int l = 0;
    /* number of bound state to be calculed */
    int n = 5;
    double dE = 0.005;
    double E;
    // double xx = h;

    double E_bound[n];
    
    /* initialization */
    Param_delta pd;
    Param_F pf;
    pd.L = L;
    pd.h = h;
    pd.n = 0; // start with the ground state
    pd.l = l;
    
    pf.l = l;
    pf.xi = 1;
    E = sqrt( (double)l*(l+1)) + dE; // minimum value of energy (?)
    printf("\nE: %lf\n\n", E);
    pf.E = E;

    FILE *file; //*fileFho3d;
    file = fopen("delta.csv","w");
    // do{
    //     fprint_double(file,E);
    //     double delta = run_for_delta3d(E,&pd);
    //     fprint_double_newline(file,delta);
    //     E += dE;
    // } while(E<10.0);

    // fileFho3d = fopen("Fho3d.csv","w");
    // do{
    //     fprint_double(fileFho3d,xx);
    //     fprint_double_newline(fileFho3d,F_ho3d(xx,&pf));
    //     xx += h;
    // } while(x<L);

     double delta_min, delta, delta_new; 

     delta_min = run_for_delta3d(E,&pd);

 do{

    delta = run_for_delta3d(E,&pd);     
    delta_new = run_for_delta3d(E+dE,&pd);
    if(delta_min * delta_new < 0.0){
        E_bound[pd.n] = bisection(run_for_delta3d,E-dE,E+dE,&pd);
        pd.n++;
        delta_min = delta_new;
        }

        E += dE;

    }while(pd.n < n);


    for(int i=0; i<n; i++){
        printf("E_%d%d = %4.3e\n",i,l,E_bound[i]);
    }
    

    
    /* eigenfunctions */
    int dim = L / h;

    double x[dim], psi[dim];

    printf("\nWhich eigenfuncion do you want to plot? (choose a number between 0 and 4) \nn = ");
    int nb; // eigenfunction to be plotted
    scanf("%d",&nb);
    assert(nb < n);


    x[0] = h;
    x[1] = x[0] + h;
    psi[0] = pow(x[0],l+1);
    psi[1] = pow(x[1],l+1);

    initialize_position(x,h,dim);

    pf.E = E_bound[nb];
    pf.l = l;
    solve_numerov(x,psi,dim,h,F_ho3d,&pf);
    normalize(psi,h,dim);

    /* save data */
    FILE *file_psi;
    file_psi = fopen("psi3d.csv","w");

    for(int i=0; i<dim; i++){
        fprint_double(file_psi,x[i]);
        fprint_double(file_psi,psi[i]);
        fprint_double_newline(file_psi,x[i]*x[i]/2.0);  // to save also the potential, but for the plot we need to increase the normalization of psi
    }
    
    
    fclose(file_psi);


}

