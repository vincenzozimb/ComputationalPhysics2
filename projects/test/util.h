#ifndef __UTIL_H__
#define __UTIL_H__

#define EPS 1e-4

/* typedef */
typedef struct Param{
    double L, h, E, R, rho;
    int l;
}Param;


/* functions */
double potential(double r, double rho, double R);
void print_potential(double L, double h, double rho, double R, int l);
double F_pot(double r, void *p);
void initialize_position(double x[], double h, int dim);
void solve_numerov(double x[], double psi[], int dim, double h, double F(double, void *), void *p);
void normalize(double psi[], double dx, int dim);
double shooting(double E, void *p);
double run_for_delta(double E, void *p);
void debug_print(double x[], double v[], int dim);

#endif
