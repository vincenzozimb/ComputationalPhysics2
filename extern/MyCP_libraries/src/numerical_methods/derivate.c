#include <assert.h>
#include <math.h>
#include <numerical_methods/derivate.h>

double derivate(double f(double,void*), double x0, double dx, void *p){
    return ( f(x0+dx,p) - f(x0-dx,p) ) / (2*dx);
}

double derivate_best(double f(double,void*), double x0, double dx, void *p){
    return ( -f(x0+2*dx,p) + 8*f(x0+dx,p) -8*f(x0-dx,p) + f(x0-2*dx,p) ) / (12*dx);
}