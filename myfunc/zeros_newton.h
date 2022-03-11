/*-------------------------------------------------------------------------*/
/**
  @file     zeros_newton_h
  @author   V. Zimbardo
  @date     Apr 2021
  @brief    C routine for finding zeros of a function using Newton-Raphson algorithm.               
*/
/*--------------------------------------------------------------------------*/


#ifndef __ZEROS_H__
#define __ZEROS_H__

#ifndef EPS
#define EPS 1e-6
#endif

/**
 * @brief Find the zero of function f starting from "guess" x_0
 * 
 * @param f function whose zero is to be calculated
 * @param x_0 initial guess 
 * @param param pointer to the struct of parameters of function f
 * 
 * @return The zero of function f
 */
double zero_newton (double f (double, void *), double x_0, void *param);

#endif