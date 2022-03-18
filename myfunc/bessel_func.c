#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "bessel_func.h"

double j_minus(double x){
    assert(x >= EPS);
    return cos(x) / x;
}

double j_zero(double x){
    assert(x >= EPS);
    return sin(x) / x;
}

double n_minus(double x){
    assert(x >= EPS);
    return j_zero(x);
}

double n_zero(double x){
    assert(x >= EPS);
    return -j_minus(x);
}

double calculate_j(double x, int order){

    assert(x >= EPS);

    double prec = j_minus(x);
    double curr = j_zero(x);

    for(int l=0; l<order; l++){
        curr = (2.0*l + 1) / x * curr - prec;
    }

    return curr;

}

double calculate_n(double x, int order){

    assert(x >= EPS);

    double prec = n_minus(x);
    double curr = n_zero(x);

    for(int l=0; l<order; l++){
        curr = (2.0*l + 1) / x * curr - prec;
    }

    return curr;

}

double recursive_bessel(double x, int l, double prec, double curr){
    assert(x >= EPS);

    return (2.0*l + 1) / x * curr - prec;
}

// double j_plus(double x){
//     assert(x >= EPS);
//     return sin(x) / (x * x) - cos(x) / x;
// }

// double n_plus(double x){
//     assert(x >= EPS);
//     return -cos(x) / (x * x) - sin(x) / x;
// }