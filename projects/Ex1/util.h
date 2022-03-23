#ifndef __UTIL_H__
#define __UTIL_H__

/* typedef */
typedef struct Param{
    double xi;
    double E;
    int l;
    double h;
    double L;
    double x0;
    int n;
}Param;


typedef struct FuncBessel{
    double j_prec, j_curr, n_prec, n_curr;
}FuncBessel;



/* functons */
double lj(double r);
double F(double r, void *par);
double Fho(double r, void *par);
double delta(double E, void *p);
void initialize_r(double r[], double dr, int dim);
void solve_numerov(double r[], double u[], int dim, double dr, double F(double, void*), void *p);
double phase_shift(double r1, double r2, double u1, double u2, double k, int l, FuncBessel *f1, FuncBessel *f2);

#endif
