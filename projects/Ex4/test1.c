#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_linalg.h>

#include "print_routines.h"
#define DIMx 1000
#define DIMe 100
#define N 8

/* ======================= FUNCTION HEADERS ======================= */

gsl_matrix *invert_a_matrix(gsl_matrix *matrix);
void print_mat_contents(gsl_matrix *matrix);
void randomize_mat_contents(gsl_matrix *matrix);
void identity_mat_contents(gsl_matrix *matrix);
void product_mat(gsl_matrix *matrix);
void initialize_mat_contents(gsl_matrix *matrix);

static size_t size = 2;

/* ======================= MAIN ======================= */
int main(){

gsl_matrix *mat = gsl_matrix_alloc(size, size);

for (size_t ii = 0; ii < size; ii++)
{
    for (size_t jj = 0; jj < size; jj++)
    {
        if (ii == jj)
        {
            gsl_matrix_set(mat, ii, jj, 1.0);
        } else {
            gsl_matrix_set(mat, ii, jj, 0.0);
        }
    } 
}

    printf("Mat2 : \n");
    print_mat_contents(mat);
    printf("\n");

gsl_matrix *inverse2 = gsl_matrix_alloc(size, size);// = invert_a_matrix(mat2);
    // printf("Inverted mat2 : \n");
    // print_mat_contents(inverse2);
    //         printf("mat[1,1] = %lf",gsl_matrix_get(mat2, 1, 1));
gsl_permutation *p = gsl_permutation_calloc(size);
int s;
    gsl_linalg_LU_decomp(mat, p, &s);
    print_mat_contents(mat);

// gsl_matrix *mat = gsl_matrix_alloc(size,size);
//     randomize_mat_contents(mat);
//     printf("Random : \n");
//     print_mat_contents(mat);
//     printf("\n");
// gsl_matrix *inverse = invert_a_matrix(mat);
//     printf("Inverted matrix:\n");
//     print_mat_contents(inverse);
gsl_matrix *prod = gsl_matrix_alloc(size, size);
printf("\n");
double val = 0.0;
    for(size_t i = 0; i < size; i++){
        for (size_t j = 0; j < size; j++){
            val = 0.0;
            for (size_t k = 0 ; k < size; k++){
                // printf("mat[%zu,%zu] = %lf", i, k ,gsl_matrix_get(mat2, i, k));
                // printf("\n");
                // printf("inverse[%zu,%zu] = %lf", k , j ,gsl_matrix_get(inverse2, k, j) );
                // printf("\n");
                // printf("mat[1,1] = %lf",gsl_matrix_get(mat2, 1, 1));

                val +=  gsl_matrix_get(mat, i, k) * gsl_matrix_get(inverse2, k, j);
                }
                printf("\n");
                gsl_matrix_set(prod, i, j , val);
        }
    }


    //gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, mat2, inverse2, 0.0, prod); 
    printf("\n");
    printf("\n");
    printf("\n");
    print_mat_contents(prod);

}
/* ======================= FUNCTION BODIES ======================= */

gsl_matrix *
invert_a_matrix(gsl_matrix *matrix)
{
    gsl_permutation *p = gsl_permutation_alloc(size);
    int s;

    // Compute the LU decomposition of this matrix
    gsl_linalg_LU_decomp(matrix, p, &s);

    // Compute the  inverse of the LU decomposition
    gsl_matrix *inv = gsl_matrix_alloc(size, size);
    gsl_linalg_LU_invert(matrix, p, inv);

    gsl_permutation_free(p);

    return inv;
}

void
print_mat_contents(gsl_matrix *matrix)
{
    size_t i, j;
    double element;

    for (i = 0; i < size; ++i) {
        for (j = 0; j < size; ++j) {
            element = gsl_matrix_get(matrix, i, j);
            printf("%f ", element);
        }
        printf("\n");
    }
}


void
randomize_mat_contents(gsl_matrix *matrix)
{
    size_t i, j;
    double random_value;
    double range = 1.0 * RAND_MAX;

    for (i = 0; i < size; ++i) {
        for (j = 0; j < size; ++j) {

            // generate a random value
            random_value = rand() / range;
    
            // set entry at i, j to random_value
            gsl_matrix_set(matrix, i, j, random_value);

        }
    }
}


void 
identity_mat_contents(gsl_matrix *matrix)
{
    size_t i,j;
    for (i = 0; i < size; i++) {
        for(j = 0; j < size; j++){

            if( i == j){
                gsl_matrix_set(matrix, i, j, 1.0);
            } else {
                gsl_matrix_set(matrix, i, j, 0.0);
            } 
        }

    }

}
