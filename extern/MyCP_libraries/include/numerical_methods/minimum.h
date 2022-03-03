/*-------------------------------------------------------------------------*/
/**
  @file     minimum_h
  @author   V. Zimbardo
  @date     Apr 2021
  @brief    C routine for calculate the minimum of a function.               
*/
/*--------------------------------------------------------------------------*/


#ifndef __MINIMUM_H__
#define __MINIMUM_H__

#ifndef EPS
#define EPS 1e-6
#endif

/**
 * @brief Find minimum of function f between given.
 *
 * Calculate the point where the function f has minimum, restricting the rearch 
 * between the given bounds.
 *
 * @param f function to find minimum of.
 * @param low lower bound for the search of the minimum.
 * @param high higher bound for the search of the minimum.
 * @param p pointer to the parameters struct of function f.
 *
 * @return double
 * The point where the function f has minimum.
 */
double minimum(double f(double, void *), double low, double high, void *p);

#endif