#include <linear_algebra/array_routines.h>

#include <assert.h>
#include <cblas.h>
#include <complex.h>
#include <math.h>
#include <stdio.h>


void mat_copy(int nrows, int ncols, const double *src, double *dest) {
	assert(nrows >= 0 && ncols >= 0);
	for (int i = 0; i < nrows; i++) {
		for (int j = 0; j < ncols; j++) {
			*(dest + i * ncols + j) = *(src + i * ncols + j);
		}
	}
}

void mat_copy_complex(int nrows, int ncols, _Complex double (*src), _Complex double (*dest)) {
	assert(nrows >= 0 && ncols >= 0);
	for (int i = 0; i < nrows; i++) {
		for (int j = 0; j < ncols; j++) {
			*(dest + i * ncols + j) = *(src + i * ncols + j);
		}
	}
}