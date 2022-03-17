/*-------------------------------------------------------------------------*/
/**
  @file     bessel_func_h
  @date     March 2022
  @brief    C-routine to calculate recursively the Bessel functions.               
*/
/*--------------------------------------------------------------------------*/

#ifndef __BESSEL_FUNC_H__
#define __BESSEL_FUNC_H__

#ifndef EPS
#define EPS 1e-6
#endif

double j_minus(double x);
double j_zero(double x);
double n_minus(double x);
double n_zero(double x);


/**
 * @brief Evaluate the spherical Bessel function of the wanted order  
 * 
 * @param x value of the independent variable
 * @param order order of the spherical Bessel function that you want to calculate
 *
 * @return double
 * The value of the spherical Bessel function of the chosen order at point x. 
 */
double calculate_j(double x, int order);

/**
 * @brief Evaluate the spherical Neumann function of the wanted order  
 * 
 * @param x value of the independent variable
 * @param order order of the spherical Neumann function that you want to calculate
 *
 * @return double
 * The value of the spherical Neumann function of the chosen order at point x. 
 */
double calculate_n(double x, int order);

#endif
