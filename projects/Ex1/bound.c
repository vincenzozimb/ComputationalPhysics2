#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "util.h"
#include "bisection.h"
#include "numerov.h"
#include "print_routines.h"

int main(){

/* ------ Parameters ------ */
double L = 6.0;
double h = 0.01;
double dE = 0.01;
double E = 2.0 * dE;
int n = 5;
int N = (int) L/h;
printf("N : %d\n", N);

/* ----- Initialization ----- */
Param par;

par.L = L;
par.h = h;
par.n = 0;
par.x0 = sqrt(2.0 * (E - dE));

double delta_prec = delta(E - dE, &par);
double delta_curr, delta_succ;
double E_bound;
FILE *pd;
pd = fopen("delta.csv", "w");
do {

 delta_curr = delta(E , &par);
 fprint_double(pd, E);
 fprint_double_newline(pd, delta_curr);
 delta_succ = delta(E + dE, &par);

if (delta_prec * delta_curr < 0.0)
{
    assert(delta(E-dE, &par) * delta(E, &par) <0.0);
    if ((delta_prec <= delta_curr && delta_curr <= delta_succ) || (delta_prec >= delta_curr && delta_curr >= delta_succ) )
    {

        E_bound = bisection( delta, E - dE, E , &par );
        printf("E_%d : %lf\n", par.n, E_bound);
        par.n++;
        delta_prec = delta_curr;
    }
    
}
E += dE;

} while(par.n <= n);

fclose(pd);

}