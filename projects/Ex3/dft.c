#include <stdio.h>
#include <stdlib.h>
#include <math.h>


int main(){

    /* system parameters */
    int N = 20;
    double rs = 3.93; // 3.93 for Na and 4.86 for K

    double R = rs * sqrt((double)N);
    double rho = 3.0 / (4.0 * M_PI * rs * rs * rs);


    /* code parameters */
    double L = 3.0 * R; // radial interval goes from 0 to L
    double h = 1e-3;

    /* number of bound state to be calculed */
    int n = 10;
    int l_max = 3;

    double E = 0.1;
    double dE = 0.05;

    double E_bound[n];




}