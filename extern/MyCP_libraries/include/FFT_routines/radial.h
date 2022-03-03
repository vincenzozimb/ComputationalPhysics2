/*-------------------------------------------------------------------------*/
/**
  @file     radial_h
  @author   V. Zimbardo
  @date     May 2021
  @brief    C routines for the radial fourier transform of a three dimensional
            radial function and its inverse. The radial fourier transform is related
            to the Hankel transform, so this code can excetute also the Hankel transform
            of order nu = 1/2. The H_{1/2} transform and its inverse are nothing more 
            than the sine and the inverse sine transforms.            
*/
/*--------------------------------------------------------------------------*/

#ifndef __RADIAL_H__
#define __RADIAL_H__

/**
 * @brief  Take a radial function and inizialize the radial fourier transform and the vector
 *         contening the magnitude of frequencies.
 * 
 * @param f the radial function
 * @param L the direct space length
 * @param N the number of points. It must be a power of 2
 * @param k the N-vector you want to inizialize with the magnitude of frequencies
 * @param F the N-vector you want to inizialize with the radial transform 
 * 
 * @return void
 */
void radial_FFT(double f[], double L, int N, double k[], double F[]);

/**
 * @brief  Take a radial fourier transform and inizialize the radial fourier inverse 
 *         transform and the vector contening the radial coordinate.
 * 
 * @param F the radial fourier transform
 * @param L the direct space length
 * @param N the number of points. It must be a power of 2
 * @param r the N-vector you want to inizialize with the radial coordinates
 * @param f the N-vector you want to inizialize with the radial inverse transform 
 * 
 * @return void
 */
void radial_IFFT(double F[], double L, int N, double r[], double f[]);

#endif