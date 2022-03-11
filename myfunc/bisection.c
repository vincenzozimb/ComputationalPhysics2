#include <assert.h>
#include <math.h>
#include "bisection.h"

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