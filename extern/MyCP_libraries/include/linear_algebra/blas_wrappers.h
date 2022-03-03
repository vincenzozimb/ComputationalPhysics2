#ifndef __BLAS_WRAPPERS_H__
#define __BLAS_WRAPPERS_H__

/**
 * @brief Copy the vector src into the vector dest
 * 
 * @param N common dimension of src and dest
 * @param src vector to be copied
 * @param dest destination vector
 * 
 * @return void
 */
void vec_copy (int N, double src[], double dest[]);

/**
 * @brief Copy the complex double vector src into the complex double vector dest
 * 
 * @param N common dimension of src and dest
 * @param src vector to be copied
 * @param dest destination vector
 * 
 * @return void
 */
void vec_complex_copy(int N, _Complex double src[], _Complex double dest[]);

/**
 * @brief Execute matrix vector multiplication
 * 
 * @param nrows number of rows of the matrix
 * @param ncols numbers of columns of the matrix
 * @param M pointer to the matrix
 * @param v pointer to the vector
 * @param res pointer to the resulting vector, res = M v
 * 
 * @return void
 */
void mat_vec_mul(int nrows, int ncols, double *M, double *v, double *res);


/**
 * @brief Execute matrix multiplication
 * 
 * @param nrows1 number of rows of the first matrix
 * @param ncols2 number of columns of the second matrix
 * @param ncols1_nrows2 number of columns of the first matrix and rows of the second one
 * @param M1 pointer to the first matrix
 * @param M2 pointer to the second matrix
 * @param res pointer to the resulting matrix, res = M1 M2
 * 
 * @return void
 */
void mat_mat_mul(int nrows1, int ncols2, int ncols1_nrows2, double *M1, double *M2, double *res);


#endif