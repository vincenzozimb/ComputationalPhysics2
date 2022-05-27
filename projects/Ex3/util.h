#ifndef __UTIL_H__
#define __UTIL_H__


/* typedef */
typedef struct ParamPot{
    double R, rho;
}ParamPot;

/* functions */
void solve_radialSE_diagonalize(int N, double r[], double v[], double E[], double psi[][N], double h, int dim);
double potential(double r, void *p);
void normalize(int N, double psi[][N], double h, int dim);
void normalize_single(double psi[], double dx, int dim);
void print_potential(double r[], double v[], int dim);
void print_wf(int N, double r[], double psi[][N], int dim, double h);
void fill_position(double r[], double h, int dim);
void fill_potential(double r[], double v[], int l, int dim, void *p);
void density(int N, double n[], double psi[][N], int dim);


#endif