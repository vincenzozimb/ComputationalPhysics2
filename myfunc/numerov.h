/*-------------------------------------------------------------------------*/
/**
  @file     numerov_h
  @author   V. Zimbardo
  @date     Apr 2021
  @brief    C function for Numerov algorithm applied to Schroedinger equation.
  In order to use this routines you have to define the appropriate F function 
  for the specific potential of interest. The algorithm can be applied to ODE
  in the form:
                  psi''(x) = F(x,p)psi(x)
  Where x is the independent variable and p is a pointer to a struct containing
  the parameters of the problem, such as the adimensional parameter
  xi = hbar^2 / 2*m*V_0*x_0^2  and the energy E in reduced units, as well as
  other possible parameter of the potential function.               
*/
/*--------------------------------------------------------------------------*/

#ifndef __NUMEROV_H__
#define __NUMEROV_H__

#include <complex.h>

/* -------------------- */ /* FUNCTION */ /*------------------------ */

/**
 * @brief Step of the Numerov algorithm to solve Schroedinger equation.
 * 
 * Performs the discrete step of the Numerov algorithm to solve psi'' = F * psi.
 * You have to know dx, F(x-dx), F(x), F(x+dx), psi(x-dx), psi(x) to obtain psi(x+dx). 
 * 
 * @param x current value of the independent variable
 * @param dx increment in the x variable
 * @param psi_curr current value of psi(x)
 * @param psi_prec value of psi(x - dx)
 * @param F function F in the equation psi''(x) = F(x)*psi(x)
 * @param p pointer to the struct of parameters
 * 
 * @return complex double
 * The value of psi(x+dx) 
 */
double complex numerov_step(double x, double dx,
                            complex double psi_curr, complex double psi_prec,
                            double F(double, void *), void *p);


/**
 * @brief Execute the Numerov algorithm and fill a vector with the values of the wave function.
 * 
 * @param x Array to store the values of the independent variable. The first two elements should be already initialized.
 * @param psi Array to store the values of the wave function. The first two elements should be already initialized.
 * @param dim Dimension of the array x and psi.
 * @param F Function F in the equation psi''(x) = F(x)*psi(x).
 * @param p Pointer to the struct with the parameters of F.
 * 
 * @return
 */
void execute_numerov_complex(double x[], double dx, complex double psi[], int dim,
                             double F(double, void *), void *p);

#endif