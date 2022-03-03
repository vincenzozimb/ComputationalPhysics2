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

#endif