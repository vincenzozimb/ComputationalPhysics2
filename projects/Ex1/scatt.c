#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "util.h"
#include "print_routines.h"
#include "bessel_func.h"

int main(){

	/* system parameters */
    double m_H = 1.00784; // u
    double m_Kr = 83.798; // u
    double mu = (m_H * m_Kr) / (m_H + m_Kr); // u
    double eps = 5.9; // meV (energy scale)
    double s = 3.18; // angstrom (length scale)


    /* conversion, do not edit */
	double amu = 931.49410242e9; // meV/c^2
    double s_ev = (s / 1.97327) * 1e-6; // ħ*c/meV
    mu = mu * amu; // meV/c^2


    /* adimensional parameter xi */
    double xi = 1.0 / (2.0 * mu * s_ev * s_ev * eps);


    /* energy in [E_start,E_end] meV */
	double E_start = 0.1 / eps;
	double E_end = 5.0 / eps;
	double dE = 0.01 / eps;


	/* other parameters */
	int l_max = 10;
	double L = 10.0;
	double dr = 0.01;    

    int dim = (int)(L / dr);
    double u[dim];
    double r[dim];

    Param par;
    par.xi = xi;


    /* initial condition */
    r[0] = 0.5;
    r[1] = r[0] + dr;
    u[0] = exp(-sqrt(4.0 / (xi*25.0) ) * pow(r[0], -5));
    u[1] = exp(-sqrt(4.0 / (xi*25.0) ) * pow(r[1], -5));


    /* two points at the end */
    int idx1 = dim-8;
    int idx2 = dim-15;


    /* initialize position */
    initialize_r(r,dr,dim);

    
    /* solving and saving */
    double k, sigma, delta;
    
    FILE* f_sigma;
    f_sigma = fopen("sigma_tot.csv","w");

    double E = E_start;
    while( E < E_end ){
        
        par.E = E;
        k = sqrt(E / xi);
        
        FuncBessel f1 = {j_minus(k*r[idx1]), j_zero(k*r[idx1]), n_minus(k*r[idx1]), n_zero(k*r[idx1])};
        FuncBessel f2 = {j_minus(k*r[idx2]), j_zero(k*r[idx2]), n_minus(k*r[idx2]), n_zero(k*r[idx2])};

        sigma = 0.0;
        for(int l=0; l<l_max; l++){

        par.l = l;
        solve_numerov(r,u,dim,dr,F,&par);
        delta = phase_shift(r[idx1],r[idx2],u[idx1],u[idx2],k,l,&f1,&f2);
        sigma += (2*l+1) * pow(sin(delta),2);
        
        }
    
        sigma *= 4*M_PI/(k*k);

        fprint_double(f_sigma,E*eps);
        fprint_double_newline(f_sigma,sigma*s*s);

        E += dE;
    }

    fclose(f_sigma);





    
  



}