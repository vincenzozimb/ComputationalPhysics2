#ifndef __UTIL_H__
#define __UTIL_H__

#define PRINT 0

/* typedef */
typedef struct Param_F{
    double xi, E;
    int l;
}Param_F;

typedef struct Param_delta{
    double L,x0,h;
    int n, l;
}Param_delta;

typedef struct FuncBessel{
    double j_prec, j_curr, n_prec, n_curr;
}FuncBessel;


/* functons */

// GENERIC
void initialize_position(double x[], double h, int dim);
void solve_numerov(double x[], double psi[], int dim, double h, double F(double, void*), void *p);

// SCATTERING 
double lj(double r);
double F_lj(double r, void *par);
double phase_shift(double r1, double r2, double u1, double u2, double k, int l, FuncBessel *f1, FuncBessel *f2);

// BOUND STATES
double F_ho(double x, void *p);
double F_ho3d(double x, void *p);
void normalize(double psi[], double dx, int dim);
void initial_condition_ho(double x[], double psi[], double L, double h, int n);
void initial_condition_ho3d(double xL[], double xR[], double psiL[], double psiR[], double L, double h, int l);
double run_for_delta(double E, void *p);
double run_for_delta3d(double E, void *p);

#endif
