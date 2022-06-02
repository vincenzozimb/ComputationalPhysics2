#ifndef __UTIL_H__
#define __UTIL_H__

#define EPS 1e-3

/* global variables */
extern int N;
extern double rs, R, rho;

/* typedef */
typedef struct ParamPot{
    double R, rho;
}ParamPot;

typedef struct ParamEc{
    double p, A, alpha1, beta1, beta2, beta3, beta4;
}ParamEc;


/* functions */
void solve_radialSE_diagonalize(int N, double r[], double v[], double E[], double psi[][N], double h, int dim);
double potential(double r, void *p);
void normalize(int N, double psi[][N], double h, int dim);
void normalize_single(double psi[], double dx, int dim, double multiplier);
void print_func(double r[], double v[], int dim, char name[25]);
void print_wf(int N, double r[], double psi[][N], int dim, double h, char name[6]);
void fill_zero(double v[], int dim);
void fill_position(double r[], double h, int dim);
void fill_potential(double r[], double v[], int l, int dim, void *p);
void add_density(int N, double r[], double n[], double psi[][N], int dim, int l);
void add_energy(int *cnt, double E[], double eps[], int Nb);
void exchange_pot(double ex[], double n[], int dim);
void add_correlation_pot(double ec[], int dim);
void add_coulomb_pot(double vh[], double r[], double n[], double h, int dim);
void partial_pot(double pot[], double r[], int dim, int l, void *p);
void change_pot(double pot[], double n[], double r[], double h, int dim);
void copy_vec(double copy[], double paste[], int dim);
double euclidean_distance(double v1[], double v2[], int dim);

#endif