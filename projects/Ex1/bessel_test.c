#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_sf_bessel.h>

#include "bessel_func.h"
#include "print_routines.h"


/* functions */
void save_gsl(int dim, int l, double dx, FILE *outfile);
void save_my(int dim, int l, double dx, FILE *outfile);

/* ------------------------------------------------------------------------------------ */
int main()
{

    // TESTING THE ROUTINE "BESSEL_FUNC" //

    /* parameters */
    double L = 15.0;
    double dx = 1e-3;

    int dim = (int)(L / dx);
    int l = 1;

    /* initialize and save */
    
    FILE *file;
    file = fopen("data_gsl.csv","w");
    save_gsl(dim,l,dx,file);
    fclose(file);

    file = fopen("data.csv","w");
    save_my(dim,l,dx,file);
    fclose(file);
    




}
/* ------------------------------------------------------------------------------------ */

/* functions */
void save_gsl(int dim, int l, double dx, FILE *outfile){

    double x,f;
    x = 1.8;

    for(int i=0; i<dim; i++){
        x += dx;
        fprint_double(outfile,x);
        for(int j=0; j<l; j++){
            f = gsl_sf_bessel_jl(j,x);
            fprint_double(outfile,f);
            f = gsl_sf_bessel_yl(j,x);
            fprint_double(outfile,f);
        }
        fprintf(outfile,"\n");
    }

}

void save_my(int dim, int l, double dx, FILE *outfile){

    double x,f;
    x = 1.8;

    for(int i=0; i<dim; i++){
        x += dx;
        fprint_double(outfile,x);
        for(int j=0; j<l; j++){
            f = calculate_j(x,j);
            fprint_double(outfile,f);
            f = calculate_n(x,j);
            fprint_double(outfile,f);
        }
        fprintf(outfile,"\n");
    }

}