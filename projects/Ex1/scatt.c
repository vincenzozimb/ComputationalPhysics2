#include  <stdio.h>
#include <stdlib.h>

#include "util.h"

int main(){

	/* system parameters */
    double m_H = 1.00784; // u
    double m_Kr = 83.798; // u
    double mu = (m_H * m_Kr) / (m_H + m_Kr); // u
    double eps = 5.9; // meV (energy scale)
    double s = 3.18; // angstrom


    /* conversion, do not edit */
	double u = 931.49410242e9; // meV/c^2
    double s_ev = (s / 1.97327) * 1e-6; // ħ*c/meV (length scale)
    mu = mu * u; // meV/c^2 (mass scale)  


    /* adimensional parameter xi */
    double xi = 1.0 / (2.0 * mu * s_ev * s_ev * eps);


    /* energy in [E_start,E_end] meV */
	double E_start = 0.5 / eps;
	double E_end = 5.0 / eps;
	double dE = 0.01 / eps;


	/* other parameters */
	int l_max = 25;
	double L = 10.0;
	double dr = 0.001;    




}