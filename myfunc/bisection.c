#include <assert.h>
#include <math.h>
#include "bisection.h"
#include "print_routines.h"

double bisection(double f(double,void*),double low, double high, void *p){
    
    assert(low < high);
	assert(isfinite(low) && isfinite(high));
	assert(f(low, p) * f(high, p) <= 0.0);

	if ( fabs(f(low, p)) <= EPS ) {
		return low;
	}
    if ( fabs(f(high, p)) <= EPS ) {
		return high;
	}

	while (fabs(high - low) / (fabs(high + low)) > EPS) {
		double mid = (high + low) / 2.0;
		if ( fabs(f(mid, p)) <= EPS) {
			return mid;
		}
		if (f(low, p) * f(mid, p) < 0.0) {
			high = mid;
		} else {
			low = mid;
		}
	}
	return low;
}

void multiple_zeros(double low, double high, double h, double f(double, void*), double zeros[], int N, void *p){

	// ADD SOMETHING FOR THE CONTROL OF THE INPUT

	int n = 0;
	double var = low;
	double f_min = f(var,p);
	double f_next;

	do{

		f_next = f(var+h,p);

		if(f_min * f_next < 0.0){
			zeros[n] = bisection(f,var-h,var+h,p);
			n++;
			f_min = f_next;
		}
		var += h;

	}while( n < N );

}