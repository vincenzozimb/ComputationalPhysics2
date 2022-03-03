/*-------------------------------------------------------------------------*/
/**
  @file     direct_wrappers_h
  @author   V. Zimbardo
  @date     May 2021
  @brief    C routine for interfacing the result of the FFT, in particular
            the gsl routine gsl_fft_complex_radix2_forward.           
*/
/*--------------------------------------------------------------------------*/


#ifndef __PRINT_ROUTINES_H__
#define __PRINT_ROUTINES_H__

#include <stdio.h>
#include <complex.h>

/*============================= DA TESTARE =================================*/

/**
 * @brief Take the vector f, the direct space length L and the numer of element of the array
 *        f (it must be a power of two). Inizialize the frequency in the vector k and
 *        the fourier transform in the vector F, in a ordered way, namely from -\piN/L to \piN/L. 
 * 
 * @param f the vector you want to transform
 * @param L direct space length
 * @param N number of point in the array (it must be a power of 2)
 * @param k empty N-vector to inizialize with the frequency
 * @param F empty N-vector to inizialize with the fourier transform
 * 
 * @return void
 */
void fourier_transform_ord(const complex double f[], double L, int N, double k[], complex double F[]);

#endif