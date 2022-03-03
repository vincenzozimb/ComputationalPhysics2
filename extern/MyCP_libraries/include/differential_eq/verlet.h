/*-------------------------------------------------------------------------*/
/**
  @file     verlet_h
  @author   V. Zimbardo
  @date     Apr 2021
  @brief    C functions for verlet algorithm. 

*/
/*--------------------------------------------------------------------------*/

#ifndef __VERLET_H__
#define __VERLET_H__

/**
 * @brief Verlet step for single particle in 1 dimension
 *
 * @param x pointer to the current position
 * @param v pointer to the current velocity
 * @param dt time increment
 * @param acc function that takes as input the position and a parameter struct, and returns the acceleration
 * @param p pointer to the struct with the parameters of acc
 *
 * @return void
 */
void verlet_1D(double *x, double *v, double dt, double acc(double, void *), void *p);

/**
 * @brief Verlet step for single particle in 2 dimensions
 *
 * @param x array with the position's components
 * @param v array with the velocity's components
 * @param a array with the acceleration's components
 * @param dt time increment
 * @param calculate_acc function that takes the position and acceleration arrays and a parameter struct, and modifies
 * the acceleration based on the current position
 * @param p pointer to struct of the parameters of calculate_acc()
 *
 * @return
 */
void verlet_2D(double x[], double v[], double a[], double dt,
			   void calculate_acc(double[], double[], void *), void *p);

/**
 * @brief Verlet step for single particle in 3 dimensions
 *
 * @param x array with the position's components
 * @param v array with the velocity's components
 * @param a array with the acceleration's components
 * @param dt time increment
 * @param calculate_acc function that takes the position and acceleration arrays and a parameter struct, and modifies
 * the acceleration based on the current position
 * @param p pointer to struct of the parameters of calculate_acc()
 *
 * @return void
 */
void verlet_3D(double x[], double v[], double a[], double dt,
			   void calculate_acc(double[], double[], void *), void *p);

/**
 * @brief Verlet step for a multibody sistem in 3 dimensions
 *
 * @param N number of bodies
 * @param x array N*3 with the positions components of the N particles
 * @param v array N*3 with the velocity components of the N particles
 * @param a array N*3 with the acceleration components of the N particles
 * @param dt time increment
 * @param calculate_acc function that takes the position and acceleration arrays and a parameter struct, and modifies
 * the acceleration based on the current position
 * @param p pointer to struct of the parameters of calculate_acc()
 * 
 * @return void
 */
void verlet_multibody_3D(int N, double x[][3], double v[][3], double a[][3], double dt,
		      			 void calculate_acc(double[][3], double[][3], void *), void *p);

/**
 * @brief Verlet step for a multibody sistem in 2 dimensions
 *
 * @param N number of bodies
 * @param x array N*2 with the positions components of the N particles
 * @param v array N*2 with the velocity components of the N particles
 * @param a array N*2 with the acceleration components of the N particles
 * @param dt time increment
 * @param calculate_acc function that takes the position and acceleration arrays and a parameter struct, and modifies
 * the acceleration based on the current position
 * @param p pointer to struct of the parameters of calculate_acc()
 * 
 * @return void
 */
void verlet_multibody_2D(int N, double x[][2], double v[][2], double a[][2], double dt,
			 			 void calculate_acc(double[][2], double[][2], void *), void *p);

#endif