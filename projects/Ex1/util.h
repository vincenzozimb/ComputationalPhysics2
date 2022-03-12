#ifndef __UTIL_H__
#define __UTIL_H__

/* typedef */
typedef struct Param{
    double xi;
    double E;
}Param;


/* functions headers */
double F_gauss(double x, void *p);
void solve_numerov(double x[], complex double psi[], int dim, double dx, double F (double, void *p), void *p, FILE *outfile);


#endif
