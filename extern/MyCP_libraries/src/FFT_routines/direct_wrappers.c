#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <assert.h>
#include <gsl/gsl_fft_complex.h>
#include <linear_algebra/blas_wrappers.h>

#include <FFT_routines/direct_wrappers.h>
/*============================= DA TESTARE =================================*/
void fourier_transform_ord(const complex double f[], double L, int N, double k[], complex double F[]) {
	
    /* check if N is a power of 2 */
    // write the assertion, using the logarithm

    /* copy f into a temporary vector */
    complex double prov[N];
    vec_complex_copy(N,f,prov);

    /* FFT of the vector F */ 
    gsl_fft_complex_radix2_forward((double *)F, 1, N);

    /* initialize k in [-pi*N/L,pi*N/L] */
	for (int n = 0; n < N; n++) {
		k[n] = 2 * M_PI / L * n - M_PI * N / L;
	}

    /* sorting of the transform vector and initialize F */
	for (int n = 0; n <= N / 2 - 2; n++) {
		F[n] = cabs(prov[N / 2 + 1 + n]);
	}
	for (int n = N / 2 - 2 + 1; n < N; n++) {
		F[n] = cabs(prov[n - N / 2 + 1]);
	}

}