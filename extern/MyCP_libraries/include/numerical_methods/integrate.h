/*-------------------------------------------------------------------------*/
/**
  @file     integrate_h
  @author   V. Zimbardo
  @date     Apr 2021
  @brief    C routine for integration.               
*/
/*--------------------------------------------------------------------------*/

#ifndef __INTEGRATE_H__
#define __INTEGRATE_H__

#ifndef EPS
#define EPS 1e-6
#endif


/**
 * @brief Evaluate the definite integral of the function f between the given limits
 * using the midpoint rule. 
 * 
 * @param f function to integrate
 * @param low lower integration limit
 * @param high higher integration limit
 * @param density point density
 * @param p pointer to the struct containing the parameters of f
 *
 * @return double
 * The value of the integral of f between low and high
 */
double integrate(double f(double, void *p), double low, double high,
                 int density, void *p);

#endif