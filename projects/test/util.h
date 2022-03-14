#ifndef __UTIL_H__
#define __UTIL_H__

/* typedef */
typedef struct ParamF{
    double xi, E;
}ParamF;

/* functions */
double pot_step(double x);
double F(double v, void *p);
void solve_numerov(double x[], double v[], complex double psi[], int dim, double dx, double F(double, void *), void *p);
void save_data(double x[], double v[], complex double psi[], int dim, FILE *outfile);


#endif
