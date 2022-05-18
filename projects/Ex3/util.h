#ifndef __UTIL_H__
#define __UTIL_H__

/* typedef */
typedef struct Param_delta{
    double L,x0,h;
    int n, l;
}Param_delta;


/* functions */
double potential(double r, double rho, double R);


#endif
