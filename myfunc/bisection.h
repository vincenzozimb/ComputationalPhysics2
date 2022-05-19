/*-------------------------------------------------------------------------*/
/**
  @file     bisection_h
  @author   V. Zimbardo
  @date     Apr 2021
  @brief    Find zero of a function using the bisection method.               
*/
/*--------------------------------------------------------------------------*/

#ifndef __BISECTION_H__
#define __BISECTION_H__

#ifndef EPS
#define EPS 1e-6
#endif


/**
 * @brief Evaluate the definite integral of the function f between the given limits
 * using the midpoint rule. 
 * 
 * @param f function to find zero
 * @param low lower limit for bisection algorithm
 * @param high higher limit for bisection algorithm
 * @param p pointer to the struct containing the parameters of f
 *
 * @return double
 * The zero of the function f between low and high (if present). 
 */
double bisection(double f(double,void*),double low, double high, void *p);



/**
 * @brief Calculate a the position of N zeros of the function f starting from the point low.
 * 
 * @param low lower limit of the interval
 * @param h ultraviolet cutoff
 * @param f function whose zeros you want to compute
 * @param zeros vector to be filled with the zeros of the function f
 * @param N maximum number of zeros to be found
 * @param p pointer to the struct containing the parameters of f
 */
void multiple_zeros(double low, double h, double f(double, void*), double zeros[], int N, void *p);

#endif