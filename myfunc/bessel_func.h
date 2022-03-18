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


/**
 * @brief Function for the recursive step for calculating the spherical Bessel and Neumann functions j_l(x) and n_l(x)
 * 
 * @param x value of the independent variable
 * @param l order of the special function s_l(x) (s could be either j or n)
 * @param prec value of s_l-1(x)
 * @param curr value of s_l(x)
 * 
 * @return The value of s_l+1(x)
 * 
 */
double recursive_bessel(double x, int l, double prec, double curr);

#endif
