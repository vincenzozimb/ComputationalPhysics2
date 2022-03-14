#include "integrate_trap.h"
#include <math.h>

double integrate(double f(double , void *), double a, double b, int density, void *param){

    int n = (int)ceil((b - a) * density);
    double dx = (b - a)/(n + 1);
    double ris = 0;
    double x1 , x2;

    for(int ii = 1; ii < n; ii++){
       x1 = a + dx * ii;
       x2 = x1 + dx;
       ris += f((x1 + x2) / 2.0, param);
    }

    return ris * dx;
}