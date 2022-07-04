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

void init_positions(double mat[DIM][N]){
    /* Initialize matrix by hand */
    mat[0][0] = 0.0;   mat[1][0] = 1.0; mat[2][0] = 0.0;
    mat[0][1] = 0.0;   mat[1][1] = 1.0; mat[2][1] = 1.0;
    mat[0][2] = 1.0;   mat[1][2] = 1.0; mat[2][2] = 2.0;
    mat[0][3] = 1.0;   mat[1][3] = 2.0; mat[2][3] = 1.0;
    mat[0][4] = 0.0;   mat[1][4] = 0.0; mat[2][4] = 0.0;
    mat[0][5] = 2.0;   mat[1][5] = 1.0; mat[2][5] = 0.0;
    mat[0][6] = 0.0;   mat[1][6] = 2.0; mat[2][6] = 0.0;
    mat[0][7] = 2.0;   mat[1][7] = 1.0; mat[2][7] = 2.0;

}

void create_submat_and_r(double mat[DIM][N], double mat_up[DIM][N2], double mat_down[DIM][N2], double r_up[N2], double r_down[N2]){
/* Create X_old and X_new matrices*/
    for (int jj = 0; jj < N; jj++){
        for (int ii = 0; ii < DIM; ii++){
            if (jj < 4){
                mat_up[ii][jj] = mat[ii][jj];
            } else {
                mat_down[ii][jj-4] = mat[ii][jj];
            }
        }
    }
    /* Create r*/
    for (int jj = 0; jj < N2; jj++){
        r_up[jj] = sqrt( mat_up[0][jj] * mat_up[0][jj] + mat_up[1][jj] * mat_up[1][jj] + mat_up[2][jj] * mat_up[2][jj] );
        r_down[jj] = sqrt( mat_down[0][jj] * mat_down[0][jj] + mat_down[1][jj] * mat_down[1][jj] + mat_down[2][jj] * mat_down[2][jj] );
    }
}

void print_matN(double mat[DIM][N]){
    for (int ii = 0; ii < DIM; ii++)
    {
        for (int jj = 0; jj < N; jj++)
        {
            printf("%lf\t", mat[ii][jj]);
        }
        printf("\n");
    }
    printf("\n");
    
}

void print_matN2(double mat[DIM][N2]){
        for (int ii = 0; ii < DIM; ii++)
    {
        for (int jj = 0; jj < N2; jj++)
        {
            printf("%lf\t", mat[ii][jj]);
        }
        printf("\n");
    }   
    printf("\n");
}

void print_vecN2(double vec[N2]){
    for (int ii = 0; ii < N2; ii++)
    {
        printf("%lf \t", vec[ii]);
    }
    printf("\n");
    
}

double Psi0(double eta, double r){
    /* wave function for orbital s */
    return  exp(-r * r /(2.0 * eta * eta));
}

double Psi1(double eta, double r, double x){
    /* wave function for orbitals p */
    return x * exp(-r * r /(2.0 * eta * eta));
}

void initialize_mat_contents(gsl_matrix *matrix, double eta, double r[N2], double R[DIM][N2]){
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
    int s;
    gsl_permutation *p = gsl_permutation_alloc(size);
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

int randomGenerator_int(int low, int up){
    /* Generates random integer numbers in the range [low, up] */
  return (rand() % (up - low + 1)) +low;
}

double randomGenerator_double(double low, double up){
    double range = up - low;
    double div = RAND_MAX / range;
    return low + (rand()/div);
}

void init_Slater_matrices(gsl_matrix *m1, gsl_matrix *m2, double eta, double r1[N2], double r2[N2], double R1[DIM][N2], double R2[DIM][N2]){
    initialize_mat_contents(m1, eta, r1, R1); // init. determinant up
    initialize_mat_contents(m2, eta, r2, R2); // init. determinant down
}

void gradient_of_matrix_x(gsl_matrix *matrix, double eta, double r[N2], double R[DIM][N2])
{   /* Calculates one component of the gradient of the matrix */
    double val = 0.0;

        for (size_t ii = 0; ii < size; ii++) // run over rows : 4 positions of electrons
        {
                val = - R[0][ii]/ (eta * eta) * Psi0(eta, r[ii]);
                gsl_matrix_set(matrix, ii, 0, val);
                
                val = ( 1.0 - R[0][ii] * R[0][ii] / (eta * eta) ) * Psi0(eta, r[ii]);
                gsl_matrix_set(matrix, ii, 1, val);

                val = - R[0][ii] * R[1][ii]/ (eta * eta) * Psi0(eta, r[ii]);
                gsl_matrix_set(matrix, ii, 2, val);
                            
                val = - R[0][ii] * R[2][ii]/ (eta * eta) * Psi0(eta, r[ii]);
                gsl_matrix_set(matrix, ii, 3, val);
                    
        }
}

void gradient_of_matrix_y(gsl_matrix *matrix, double eta, double r[N2], double R[DIM][N2])
{   /* Calculates one component of the gradient of the matrix */
    double val = 0.0;

        for (size_t ii = 0; ii < size; ii++) // run over rows : 4 positions of electrons
        {
                val = - R[1][ii]/ (eta * eta) * Psi0(eta, r[ii]);
                gsl_matrix_set(matrix, ii, 0, val);
                
                val = - R[0][ii] * R[1][ii]/ (eta * eta) * Psi0(eta, r[ii]);
                gsl_matrix_set(matrix, ii, 1, val);

                val = ( 1.0 - R[1][ii] * R[1][ii] / (eta * eta) ) * Psi0(eta, r[ii]);
                gsl_matrix_set(matrix, ii, 2, val);
                            
                val = - R[1][ii] * R[2][ii]/ (eta * eta) * Psi0(eta, r[ii]);
                gsl_matrix_set(matrix, ii, 3, val);
                    
        }
}

void gradient_of_matrix_z(gsl_matrix *matrix, double eta, double r[N2], double R[DIM][N2])
{   /* Calculates one component of the gradient of the matrix */
    double val = 0.0;

        for (size_t ii = 0; ii < size; ii++) // run over rows : 4 positions of electrons
        {
                val = - R[2][ii]/ (eta * eta) * Psi0(eta, r[ii]);
                gsl_matrix_set(matrix, ii, 0, val);
                
                val = - R[0][ii] * R[2][ii]/ (eta * eta) * Psi0(eta, r[ii]);
                gsl_matrix_set(matrix, ii, 1, val);

                val = - R[1][ii] * R[2][ii]/ (eta * eta) * Psi0(eta, r[ii]);
                gsl_matrix_set(matrix, ii, 2, val);
                            
                val = ( 1.0 - R[2][ii] * R[2][ii] / (eta * eta) ) * Psi0(eta, r[ii]);
                gsl_matrix_set(matrix, ii, 3, val);
                    
        }
}

void laplacian_of_matrix(gsl_matrix *matrix, double eta, double r[N2], double R[DIM][N2])
{   /* Calculated the laplacian of the matrix */
    double val;
      
        for (size_t ii = 0; ii < size; ii++)
        {   val = 0.0;
                val = Psi0(eta, r[ii]) / (eta * eta) * (-3.0 + r[ii] * r[ii] / (eta * eta) );
                gsl_matrix_set(matrix, ii, 0, val);
                val = Psi0(eta, r[ii]) * R[0][ii] / (eta * eta) * (-5.0 + r[ii] * r[ii]/ (eta * eta) );
                gsl_matrix_set(matrix, ii, 1, val);
                val = Psi0(eta, r[ii]) * R[1][ii] / (eta * eta) * (-5.0 + r[ii] * r[ii]/ (eta * eta) );
                gsl_matrix_set(matrix, ii, 2, val);
                val = Psi0(eta, r[ii]) * R[2][ii] / (eta * eta) * (-5.0 + r[ii] * r[ii]/ (eta * eta) );
                gsl_matrix_set(matrix, ii, 3, val);
        }
    
}
    
void GDtoDR(gsl_matrix *m1, gsl_matrix *m2, gsl_matrix *m3 , gsl_matrix *detInv, double GDratio[DIM][N2])
{   /* Calculates the gradient determinant-to-determinant ratio when no move is performed */
    double val = 0.0, inv = 0.0;
    double totx= 0.0, toty = 0.0, totz = 0.0;
    for (size_t ii = 0; ii < 4; ii++) // run over rows : position of electrons
        {  
            totx = toty = totz = 0.0;
            for (size_t jj = 0; jj < size; jj++) // run over column : electrons
            {   // product row by row
                //  printf("val[%zu][%zu]: %lf \t inv[%zu][%zu]: %lf\t v*i: %lf\t tot: %lf\n", ii, jj, val, ii, jj, inv, val*inv, totz);
                val =  gsl_matrix_get(m1, ii, jj);  inv = gsl_matrix_get(detInv, jj, ii );
                totx += val *inv;
                GDratio[0][ii] = totx; // x-component
                
                val =  gsl_matrix_get(m2, ii, jj);
                toty += val *inv;
                GDratio[1][ii] = toty; // y-component
                
                val = gsl_matrix_get(m3, ii, jj);
                totz += val * inv;
                GDratio[2][ii] = totz; // z-component
            }
  
    }

}

void LDtoDR(gsl_matrix *mm, gsl_matrix *detInv, double LDratio[N2])
{   /* Calculated the laplacian determinant-to determinant ratio */
    double val, inv, tot;
    for (size_t ii = 0; ii < size; ii++) // run over rows
    {   tot = 0.0;
        for (size_t jj = 0; jj < size; jj++) // run over columns
        {   val = gsl_matrix_get(mm, ii, jj);   inv = gsl_matrix_get(detInv, jj, ii);
            tot += val*inv;
        }
    LDratio[ii] = tot;
    }
    
}

double localEnergy1( double LDtoDRup[N2], double LDtoDRdown[N2]){
    /* First way to implement the calculation of the local kinetic energy */
   double val = 0.0;
    for (int ii = 0; ii < N2; ii++)
    {
             val += LDtoDRup[ii] +  LDtoDRdown[ii];
    }
    return -1.0/2.0 * val;
}

double localEnergy1_b(gsl_matrix *lapUp, gsl_matrix *lapDown, gsl_matrix *detIup, gsl_matrix *detIdown, gsl_matrix *gradUpx, gsl_matrix *gradUpy, gsl_matrix *gradUpz, gsl_matrix *gradDownx, gsl_matrix *gradDowny, gsl_matrix *gradDownz){
    double val = 0.0, lap = 0.0, detI = 0.0, val2 = 0.0;;
    for (size_t ii = 0; ii < N; ii++)
    {
        if (ii < N2)
        {
            for (int jj = 0; jj < N2; jj++)
            {
                lap = gsl_matrix_get(lapUp, ii, jj);
                detI = gsl_matrix_get(detIup, ii, jj);
                val += lap * detI;
                val2 += gsl_matrix_get(gradUpx, 0, jj) * gsl_matrix_get(gradUpx, 0, jj) + gsl_matrix_get(gradUpx, 1, jj) * gsl_matrix_get(gradUpx, 1, jj) + gsl_matrix_get(gradUpx, 2, jj) * gsl_matrix_get(gradUpx, 2, jj);
            } 
            
        } else {
            for (size_t jj = 0; jj < N2; jj++)
            {
                lap = gsl_matrix_get(lapDown, ii-4, jj);
                detI = gsl_matrix_get(detIdown, ii-4, jj);
                val += lap * detI;
                val2 += gsl_matrix_get(gradDownx, 0, jj) * gsl_matrix_get(gradDownx, 0, jj) + gsl_matrix_get(gradDownx, 1, jj) * gsl_matrix_get(gradDownx, 1, jj) + gsl_matrix_get(gradUpx, 2, jj) * gsl_matrix_get(gradDownx, 2, jj);
            }
        }
    
    }
    return -0.25 * (val - val2);
}

double localEnergy2(double LDtoDRup[N2], double LDtoDRdown[N2], double GDtoDRup[DIM][N2], double GDtoDRdown[DIM][N2])
{ /* Second way to implement the calculation of the local kinetic energy */

    double val1 = 0.0;
    double val2 = 0.0;
    for (int ii = 0; ii < N2; ii++)
    {
             val1 += LDtoDRup[ii] +  LDtoDRdown[ii];
    }

    for (int jj = 0; jj < N; jj++)
    {  
        if (jj < N2)
        {
            val2 += GDtoDRup[0][jj] * GDtoDRup[0][jj]  + GDtoDRup[1][jj] * GDtoDRup[1][jj] + GDtoDRup[2][jj] * GDtoDRup[2][jj]; 
        } else {
            val2 += GDtoDRdown[0][jj-4] * GDtoDRdown[0][jj-4]  + GDtoDRdown[1][jj-4] * GDtoDRdown[1][jj-4] + GDtoDRdown[2][jj-4] * GDtoDRdown[2][jj-4]; 
        }
        
    }
 
    return 1.0 / 4.0 * (- val1 + val2);
}

double determinant_of_matrix(gsl_matrix *matrix)
{   /* returns the Slater determinant of the input matrix */
    double detMat = 0.0;
    int s;
    gsl_permutation *p = gsl_permutation_alloc(size);
 
    // Compute the LU decomposition of this matrix (s = 1)
    gsl_linalg_LU_decomp(matrix, p, &s);
    detMat = gsl_linalg_LU_det(matrix, s);

    gsl_permutation_free(p);
    return detMat;
}

void getRatio(gsl_matrix *matrix, gsl_matrix *detInv, double ratio[N2] ){
   
    double val = 0.0; 
    double inv = 0.0;
    double tot = 0.0;
    for (size_t ii = 0; ii < size; ii++)
    {   tot = 0.0;
        for (size_t jj = 0; jj < size; jj++)
        {
            val = gsl_matrix_get(matrix, ii, jj);   inv = gsl_matrix_get(detInv, ii, jj);
            tot += val*inv;
            printf("val: %lf\t inv: %lf\t tot: %lf\n", val, inv, tot);
        }
        ratio[ii] = tot;
    }
    
}
