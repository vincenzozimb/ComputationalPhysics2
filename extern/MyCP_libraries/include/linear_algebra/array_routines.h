#ifndef __ARRAY_ROUTINES_H__
#define __ARRAY_ROUTINES_H__

#include <stdio.h>
#include <complex.h>


/**
 * @brief Copy a multimentional array.
 * 
 * Copy a multimentional array from src to dest.
 * 
 * @param nrows Number of rows of the multidimentional array. It must be >= 0.
 * @param ncols Number of columns of the multidimentional array. It must be >= 0.
 * @param src Pointer to the source array that must be copied.
 * @param dest Pointer to the destination array where src must be copied.
 * 
 * @return
 */
void mat_copy(int nrows, int ncols, const double* src, double* dest);

/**
 * @brief Copy a multimentional array of complex numbers.
 * 
 * Copy a multimentional array of comlex numbers from src to dest.
 * 
 * @param nrows Number of rows of the multidimentional array. It must be >= 0.
 * @param ncols Number of columns of the multidimentional array. It must be >= 0.
 * @param src Pointer to the source array that must be copied.
 * @param dest Pointer to the destination array where src must be copied.
 * 
 * @return
 */
void mat_copy_complex(int nrows, int ncols, _Complex double (*src), _Complex double (*dest));

#endif