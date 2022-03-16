#ifndef __UTIL_H__
#define __UTIL_H__

/* typedef */
typedef struct ParamF{
    double xi, E;
}ParamF;

/* functions */
double potential(double x);
double F_step(double x, void *p);
void solve_numerov(double x[], complex double psi[], double L, double dx, double F(double, void *), void *p);
void fill_potential(double v[], double x[], int dim);
void save_data(double x[], double v[], complex double psi[], int dim, FILE *outfile);

#endif
