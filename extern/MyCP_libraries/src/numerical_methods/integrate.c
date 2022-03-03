#include <assert.h>
#include <math.h>
#include <numerical_methods/integrate.h>

double integrate(double f(double, void *p), double low, double high,
                 int density, void *p){

	assert(density > 0);
	assert(isfinite(low) && isfinite(high));
	assert(low <= high);

	int n = 1 + (int)ceil((high - low) * density); /* number of points */
	double dx = (high - low) / n;	   /* increment in integration variable */
	double ris = 0;

	for (int i=0;i<n;i++){
		double x1 = low + i*dx;
		double x2 = x1 + dx;
		double mean = ( x1 + x2 ) / 2.0;
		ris += f(mean,p);
	}
	ris *= dx;
	
	return ris;

}