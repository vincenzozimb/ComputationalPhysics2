#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <string.h>
#include <time.h>

#include "util.h"

/* ======================= FUNCTION BODIES ======================= */

double Psi0(double eta, double r){
    /* wave function for orbital s */
    return exp(-r * r /(2.0 * eta * eta));
}

double Psi1(double eta, double r, double x){
    /* wave function for orbitals p */
    return x * exp(-r * r /(2.0 * eta * eta));
}

void initialize_mat_contents(gsl_matrix *matrix, double eta, double r[4], double R[3][4]){
    /* initialize Slater determinant matrix with single-electron wave functions */
    double val = 0.0; 
    for ( size_t jj = 0; jj < size; jj++) // run over columns: 4 electrons
    {
        for ( size_t ii = 0; ii < size; ii++) // run over rows : 4 positions of electrons
        {
            if(jj == 0){
                val = Psi0(eta, r[ii]);
                gsl_matrix_set(matrix, ii, jj, val);
            } else {
                val = Psi1(eta, r[ii], R[jj-1][ii]);
                gsl_matrix_set(matrix, ii, jj, val);
            }

        }
    
    }
    
}

gsl_matrix *invert_a_matrix(gsl_matrix *matrix)
{   /* returns the inverse of the input matrix */
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

void print_mat_contents(gsl_matrix *matrix)
{   /* prints the content of a square gsl matrix */
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

double determinant_of_matrix(gsl_matrix *matrix)
{   /* returns the Slater determinant of the input matrix */
    double detMat = 0.0;
        
    gsl_permutation *p = gsl_permutation_alloc(size);
    int s;

    // Compute the LU decomposition of this matrix (s = 1)
    gsl_linalg_LU_decomp(matrix, p, &s);
    detMat = gsl_linalg_LU_det(matrix, s);

    gsl_permutation_free(p);
    return detMat;
}

double Vext(double r, double R, double rho)
{   /* External potential */
    double a = 2.0 * M_PI * rho;
    if (r <= R)
    {
        return a * (r * r/3.0 - R * R);
    } else {
        return a * -2.0/3.0 * R * R * R/r;
    
    }
    
}

int randomGenerator(int low, int up){
    /* Generates random integer numbers in the range [low, up] */
  return (rand() % (up - low + 1)) +low;
}

void gradient_of_matrix(gsl_matrix *matrix, double eta, double r[4], double R[3][4], int component)
{   /* Calculates one component of the gradient of the matrix */
    double val = 0.0;

    for (size_t jj = 0; jj < size; jj++) // run over columns : electrons
    {
        for (size_t ii = 0; ii < size; ii++) // run over rows : 4 positions of electrons
        {
            if (jj == 0 )
            {
                val = - R[component][ii] / (eta * eta) * Psi0(eta, r[ii]);
                gsl_matrix_set(matrix, ii, jj, val);
            } else {
                val = (1.0 - (R[component][ii] * R[component][ii])  / (eta * eta)) * Psi0(eta, r[ii]);
                gsl_matrix_set(matrix, ii, jj, val);
            }
            
        }
        
    }
}
    
void laplacian_of_matrix(gsl_matrix *matrix, double eta, double r[4], double R[3][4] )
{   /* Calculated the laplacian of the matrix */
    double val;
    for (size_t jj = 0; jj < size; jj++)
    {   
        for (size_t ii = 0; ii < size; ii++)
        {   val = 0.0;
            if (jj == 0)
            {
                for (int cc = 0; cc < DIM; cc++)
                {
                    val += (- 1.0 / (eta * eta) + pow(R[cc][ii]/ (eta * eta), 2.0)) * Psi0(eta, r[ii]);
                    gsl_matrix_set(matrix, ii, jj, val);
                } 
            } else {
                for (int cc = 0; cc < DIM; cc++)
                {
                    val += ( -3.0 * R[cc][ii] / (eta * eta) + pow(R[cc][ii] , 3.0 )/ pow(eta, 4.0 ) ) * Psi0(eta, r[ii]);
                    gsl_matrix_set(matrix, ii, jj, val);
                }
    
            }
            
        }
        
    }
    
}

void GDtoDR_old(gsl_matrix *m1, gsl_matrix *m2, gsl_matrix *m3 , gsl_matrix *detInv, double GDratio[DIM][N2])
{   /* Calculates the gradient determinant-to-determinant ratio when no move is performed */
    
    for (size_t ii = 0; ii < size; ii++) // run over rows : position of electrons
        {
            for (size_t jj = 0; jj < size; jj++) // run over column : electrons
            {   // product row by row
                GDratio[0][ii] += gsl_matrix_get(m1, ii, jj) * gsl_matrix_get(detInv, ii, jj ); // x-component
                GDratio[1][ii] += gsl_matrix_get(m2, ii, jj) * gsl_matrix_get(detInv, ii, jj ); // y-component
                GDratio[2][ii] += gsl_matrix_get(m3, ii, jj) * gsl_matrix_get(detInv, ii, jj ); // z-component
            }
    }

}

void GDtoDR_new(gsl_matrix *A, gsl_matrix *m1, gsl_matrix *m2, gsl_matrix *m3 , gsl_matrix *detInv, double GDratio[DIM][N2])
{   /* Calculates the gradient determinant-to-determinant ratio when a move is performed */
    double Rsd = 0.0;
    for (size_t ii = 0; ii < size; ii++) // run over rows : position of electrons
        {
            Rsd = 0.0;
            for (size_t jj = 0; jj < size; jj++)
            {
                Rsd += gsl_matrix_get(A, ii, jj) * gsl_matrix_get(detInv, ii, jj );
            }
            for (size_t jj = 0; jj < size; jj++) // run over column : electrons
            {   // product row by row
                GDratio[0][ii] += gsl_matrix_get(m1, ii, jj) * gsl_matrix_get(detInv, ii, jj )/Rsd; // x-component
                GDratio[1][ii] += gsl_matrix_get(m2, ii, jj) * gsl_matrix_get(detInv, ii, jj )/Rsd; // y-component
                GDratio[2][ii] += gsl_matrix_get(m3, ii, jj) * gsl_matrix_get(detInv, ii, jj )/Rsd; // z-component
            }
    }

}

void LDtoDR(gsl_matrix *mm, gsl_matrix *detInv, double LDratio[N2])
{   /* Calculated teh laplacian determinant-to determinant ratio */
    
    for (size_t ii = 0; ii < size; ii++) // run over rows
    {
        for (size_t jj = 0; jj < size; jj++) // run over columns
        {
            LDratio[ii] += gsl_matrix_get(mm, ii, jj) * gsl_matrix_get(detInv, ii, jj);
        }
        
    }
    
}

void initialization_of_wf(gsl_matrix *m, gsl_matrix *gm_x, gsl_matrix *gm_y, gsl_matrix *gm_z, gsl_matrix *lm, double GDrat[DIM][N2] , double LDrat[4], double detm, double eta, double r[4], double R[3][4])
{
    gsl_matrix *inv_m;
    initialize_mat_contents(m, eta, r, R); //Initialize matrix m
    gradient_of_matrix(gm_x, eta, r, R, 0); // Gradient of m - gm_x
    gradient_of_matrix(gm_y, eta, r, R, 1); // Gradient of m - gm_y
    gradient_of_matrix(gm_z, eta, r, R, 2); // Gradient of m - gm_z
    laplacian_of_matrix(lm, eta, r, R); // Laplacian of m
    inv_m = invert_a_matrix(m); // Inverse of m
    initialize_mat_contents(m, eta, r, R); // Re-initialize m
    GDtoDR_old(gm_x, gm_y, gm_z, inv_m, GDrat); 
    LDtoDR(lm, inv_m, LDrat); // LDtoD ratio
    detm = determinant_of_matrix(m); // Slater determinant
    initialize_mat_contents(m, eta, r, R); // Re-initialize m
}

void evolution_of_wf(gsl_matrix *m, gsl_matrix *inv_m, gsl_matrix *gm_x, gsl_matrix *gm_y, gsl_matrix *gm_z, gsl_matrix *lm, double GDrat[DIM][N2] , double LDrat[4], double detm, double eta, double r[4], double R[3][4])
{
    
    initialize_mat_contents(m, eta, r, R); //Initialize matrix m
    gradient_of_matrix(gm_x, eta, r, R, 0); // Gradient of m - gm_x
    gradient_of_matrix(gm_y, eta, r, R, 1); // Gradient of m - gm_y
    gradient_of_matrix(gm_z, eta, r, R, 2); // Gradient of m - gm_z
    laplacian_of_matrix(lm, eta, r, R); // Laplacian of m
    inv_m = invert_a_matrix(m); // Inverse of m
    initialize_mat_contents(m, eta, r, R); // Re-initialize m
    GDtoDR_new(m, gm_x, gm_y, gm_z, inv_m, GDrat); 
    LDtoDR(lm, inv_m, LDrat); // LDtoD ratio
    detm = determinant_of_matrix(m); // Slater determinant
    initialize_mat_contents(m, eta, r, R);
}

double localEnergy1( double LDtoDRup[N2], double LDtoDRdown[N2])
{ /* First mode to implement the calculation of the local kinetic energy */
    double val = 0.0;
    for (int ii = 0; ii < N2; ii++)
    {
        val += LDtoDRup[ii] + LDtoDRdown[ii];
    }
    
   return -1.0 / 2.0 * val;
    
}

double localEnergy2(double LDtoDRup[N2], double LDtoDRdown[N2], double GDtoDRup[DIM][N2], double GDtoDRdown[DIM][N2])
{ /* Second way to implement the calculation of the local kinetic energy */

    double val1 = 0.0;
    double val2 = 0.0;
    for (int ii = 0; ii < N2; ii++)
    {
        val1 += LDtoDRup[ii] + LDtoDRdown[ii];

    }
    for (int jj = 0; jj < N2; jj++) // columns (1,2,3,4)
    {
        for (int ii = 0; ii < DIM; ii++) // rows (x,y,z)
        {
           val2 += pow(GDtoDRup[ii][jj] + GDtoDRdown[ii][jj], 2.0);
        }
        
    }
    
    return 1.0 / 4.0 * (- val1 + val2);
}
