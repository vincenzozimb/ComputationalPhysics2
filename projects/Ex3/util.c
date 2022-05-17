#include <stdio.h>
#include <stdlib.h>
#include <math.h>


/* functions */
double potential(double r, double rho, double R){

    double pot;
    if(r <= R){
        r * r / 3.0 - R * R;
    }else if(r > R){
        -2.0/3.0 * R * R * R / r; 
    }

    return 2.0 * M_PI * rho * pot;

}