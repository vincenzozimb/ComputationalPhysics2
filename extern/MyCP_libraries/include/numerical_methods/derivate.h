/*-------------------------------------------------------------------------*/
/**
  @file     derivate_h
  @author   V. Zimbardo
  @date     Apr 2021
  @brief    C routines for numerical derivation.               
*/
/*--------------------------------------------------------------------------*/

#ifndef __DERIVATE_H__
#define __DERIVATE_H__

/**
 * @brief Evaluate the first derivative of the function f in the point x0, using finite
 * difference method of order O(dx^2). 
 * 
 * @param f function you want to calculate the first derivative
 * @param x0 point where you want to calculate the derivative of the function f
 * @param dx accuracy for finite difference method
 * @param p pointer to the struct containing the parameters of f
 *
 * @return double
 * The value of the derivative of f in the point x0, with an accuracy of order O(dx^2)
 */
double derivate(double f(double,void*), double x0, double dx, void *p);

/**
 * @brief Evaluate the first derivative of the function f in the point x0, using finite
 * difference method of order O(dx^4). 
 * 
 * @param f function you want to calculate the first derivative
 * @param x0 point where you want to calculate the derivative of the function f
 * @param dx accuracy for finite difference method
 * @param p pointer to the struct containing the parameters of f
 *
 * @return double
 * The value of the derivative of f in the point x0, with an accuracy of order O(dx^4)
 */
double derivate_best(double f(double,void*), double x0, double dx, void *p);

#endif