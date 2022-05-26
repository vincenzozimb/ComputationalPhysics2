#ifndef __UTIL_H__
#define __UTIL_H__


/* typedef */
typedef struct ParamPot{
    double R, rho;
}ParamPot;

/* functions */
void solve_radialSE_diagonalize(int N, int l, double E[], double psi[][N], double L, int dim, double pot(double, void*), void *p, int bol);
double potential(double r, void *p);
double harmonic(double r, void *p);
void normalize(int N, double psi[][N], double h, int dim);
void normalize_single(double psi[], double dx, int dim);
void print_wf(int N, double psi[][N], int dim, double h);


#endif