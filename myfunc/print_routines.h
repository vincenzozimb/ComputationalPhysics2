/*-------------------------------------------------------------------------*/
/**
  @file     print_routines_h
  @author   V. Zimbardo
  @date     Apr 2021
  @brief    C routines for print.               
*/
/*--------------------------------------------------------------------------*/


#ifndef __PRINT_ROUTINES_H__
#define __PRINT_ROUTINES_H__

#include <stdio.h>
#include <complex.h>

/**
 * @brief Print a double to a specified file, followed by a tab.
 * 
 * @param file the file handler of type FILE for the ouput file. If you want to print to terminal, use stdout
 * @param d the double to be printed to file
 * 
 * @return void
 */
void fprint_double (FILE *file, double d);

/**
 * @brief Print a double to a specified file, followed by a newline.
 * 
 * @param file The file handler of type FILE for the ouput file. If you want to print to terminal, use stdout
 * @param d The double to be printed to file
 * 
 * @return
 */
void fprint_double_newline (FILE *file, double d);

/**
 * @brief Print a double vector to a specified file.
 * 
 * @param file the file handler of type FILE for the ouput file. If you want to print to terminal, use stdout
 * @param vec the vector of double to be printed to file
 * 
 * @return void
 */
void fprint_vec(FILE *file, double vec[], int dim);

/**
 * @brief Print multidimentional array to specified output.
 * 
 * Print an array with size nrows*ncols to the specified output.
 * 
 * @param file FILE pointer to the intended output file. To print to standard output, type stdout.
 * @param A Pointer to the multidimentional array
 * @param nrows Number of rows of the multidimentional array. It must be >= 0.
 * @param ncols Number of columns of the multidimentional array. It must be >= 0.
 * 
 * @return
 */
void fprint_mat(FILE *file, const double *A, int nrows, int ncols);

/**
 * @brief Print multidimentional array of complex numbers to specified output.
 * 
 * Print an array of complex numbers with size nrows*ncols to the specified output.
 * 
 * @param file FILE pointer to the intended output file. To print to standard output, type stdout.
 * @param A Pointer to the multidimentional array
 * @param nrows Number of rows of the multidimentional array. It must be >= 0.
 * @param ncols Number of columns of the multidimentional array. It must be >= 0.
 * 
 * @return
 */
void fprint_mat_complex(FILE *file, _Complex double (*A), int nrows, int ncols);

/**
 * @brief Print two vectors in two columns to specified output.
 * 
 * Print the two vectors v1 and v2 to the specified output.
 * 
 * @param file FILE pointer to the intended output file. To print to standard output, type stdout.
 * @param v1 pointer to the first vector
 * @param v2 pointer to the second vector
 * @param dim common dimension of the two vectors
 * 
 * @return void 
 */
void fprint_two_vec(FILE *file, double v1[], double v2[], int dim);

/**
 * @brief Print three vectors in three columns to specified output.
 * 
 * Print the three vectors v1,v2 and v3 to the specified output.
 * 
 * @param file FILE pointer to the intended output file. To print to standard output, type stdout.
 * @param v1 pointer to the first vector
 * @param v2 pointer to the second vector
 * @param v3 pointer to the third vector
 * @param dim common dimension of the three vectors
 * 
 * @return void 
 */
void fprint_three_vec(FILE *file, double v1[], double v2[], double v3[], int dim);


#endif

