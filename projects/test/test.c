#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>

#include "bisection.h"
#include "print_routines.h"

/* typedef */


/* functions */
double func(double x, void *p);
void fill_zero(double x[], int dim);
void fill_one(double x[], int dim);

/* ------------------------------------------- */
int main(){

    double L = 12.2;
    double h = 1e-2;

    // int dim = (int)(L/h);

    void *p;

    // double x;
    // FILE *file;
    // file = fopen("data.csv","w");
    // for(int i=0; i<dim; i++){
        
    //     x = (i+1) * h;
    //     fprint_double(file,x);
    //     fprint_double_newline(file, func(x,&p));

    // }
    // fclose(file);

    int n = 10;
    double zeros[n];
    fill_zero(zeros,n);

    multiple_zeros(h,L-h,h,func,zeros,n,&p);
    fprint_vec(stdout,zeros,n);


}
/* ------------------------------------------- */

/* functions */
double func(double x, void *p){
    return sin(M_PI * x) / x;
}

void fill_zero(double x[], int dim){
    for(int i=0; i<dim; i++){
        x[i] = 0.0;
    }
}

void fill_one(double x[], int dim){
    for(int i=0; i<dim; i++){
        x[i] = 1.0;
    }
}