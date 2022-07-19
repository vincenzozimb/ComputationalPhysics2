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

double Vext(double rUp[N2], double rDown[N2], double R, double rho)
{   /* External potential */
    double a = 2.0 * M_PI * rho;
    double pot = 0.0;

    for (int jj = 0; jj < N; jj++)
    {
        if (jj < N2)
        {
            if (rUp[jj] < R)
            {
                pot += a * (rUp[jj] * rUp[jj]/3.0 - R * R);
            } else {
                pot += a * -2.0/3.0 * R * R * R/rUp[jj];
            }
            
        } else {
            if (rDown[jj-4] < R)
            {
                pot += a * (rDown[jj-4] * rDown[jj-4]/3.0 - R * R);
            } else {
                pot += a * -2.0/3.0 * R * R * R/rDown[jj-4];
            }
            
        }
        
    }
    return pot;
    
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

// -----------------------------------------------------------------------------------
void init_positions(double mat[DIM][N], double delta){

    for (int ii = 0; ii < DIM; ii++)
    {
        for (int jj = 0; jj < N; jj++)
        {
            mat[ii][jj] = delta * (rand()/(double)RAND_MAX - 1.0/2.0);
        }
        
    }
    
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

double Psi0(double eta, double r){
    /* wave function for orbital s */
    return  exp(-r * r /(2.0 * eta * eta));
}

double Psi1(double eta, double r, double x){
    /* wave function for orbitals p */
    return x * exp(-r * r /(2.0 * eta * eta));
}

double orbital(double R[DIM][N2], double r[N2], double eta, int ii, int jj){
// jj -> column : orbital / electron
// ii -> row : electron position
    if (jj == 0){
        return Psi0(eta, r[ii]);
    } else {
        return Psi1(eta, r[ii], R[jj-1][ii]);
    }
    
}

void initialize_mat_contents(gsl_matrix *matrix, double eta, double r[N2], double R[DIM][N2]){
    /* initialize Slater determinant matrix with single-electron wave functions */
    for ( size_t jj = 0; jj < size; jj++) // run over columns: 4 electrons
    {
        for ( size_t ii = 0; ii < size; ii++) // run over rows : 4 positions of electrons
        {
            gsl_matrix_set(matrix, ii, jj, orbital(R, r, eta, ii, jj));

        }
    
    }
    
}

void init_Slater_matrices(gsl_matrix *m1, gsl_matrix *m2, double eta, double r1[N2], double r2[N2], double R1[DIM][N2], double R2[DIM][N2]){
    initialize_mat_contents(m1, eta, r1, R1); // init. determinant up
    initialize_mat_contents(m2, eta, r2, R2); // init. determinant down
}

int k_delta(int i, int j)
{
    if(i==j) return 1;
    else return 0;
}

double grad_orbital(double R[DIM][N2], double r[N2], double eta, int ii, int jj, int kk){
// kk -> x, y, z (1,2,3)
// jj -> orbital / electron
// ii -> electron position 
    if( jj == 0)
    {
        return -R[kk][ii]/ (eta * eta) * orbital(R, r, eta, ii, 0);
    }
    else
    {
        return orbital(R, r, eta, ii, 0) * (k_delta(jj, kk+1) - R[kk][ii] * R[jj-1][ii]/pow(eta, 2.0));
    }

}

void initialize_grad_mat(gsl_matrix* grad_mat[], double R[DIM][N2], double r[N2], double eta)
{
    for(int kk = 0; kk < DIM; kk++)
    {
        for(int ii = 0; ii < N2; ii++)
        {
            for(int jj = 0; jj < N2; jj++)
                gsl_matrix_set(grad_mat[kk], ii, jj, grad_orbital(R, r, eta, ii, jj, kk));
        }
    }
}

double laplacian_orbital(double R[DIM][N2], double r[N2] ,int ii, int jj, double eta)
{
    if(jj == 0)
    {
        return  orbital(R, r, eta, ii, 0)* (-3.0 + r[ii] * r[ii] / (eta * eta) ) / (eta * eta); 
    }
    else
    {
        return orbital(R, r, eta, ii, jj) * (-5.0 + r[ii] * r[ii]/ (eta * eta) ) / (eta * eta); 
    }
}

void initialize_lap_mat(gsl_matrix* mat, double R[DIM][N2], double r[N2], double eta)
{
    for(int ii = 0; ii < N2; ii++)
    {
        for(int jj = 0; jj < N2; jj++)
        {
            gsl_matrix_set( mat, ii, jj, laplacian_orbital(R, r, ii, jj, eta) );
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
            printf("%e \t", element);
        }
        printf("\n");
    }
}

double localEnergy1(gsl_matrix* Aup, gsl_matrix* Adown, double Rup[DIM][N2], double Rdown[DIM][N2], double r1[N2], double r2[N2], double eta){
    /* First way to implement the calculation of the local kinetic energy */
   
    gsl_matrix* lapl_Aup = gsl_matrix_alloc(N2 , N2);
    initialize_lap_mat(lapl_Aup, Rup, r1, eta);
    gsl_matrix* lapl_Adown = gsl_matrix_alloc(N2, N2);
    initialize_lap_mat(lapl_Adown, Rdown, r2, eta);

    gsl_matrix* inv_Aup = invert_a_matrix(Aup);
    gsl_matrix* inv_Adown = invert_a_matrix(Adown);

    double sum = 0;
    for(int ii = 0; ii < N2; ii++)
    {
        for(int jj = 0; jj < N2; jj++)
        {
            sum += gsl_matrix_get(lapl_Aup, ii, jj) * gsl_matrix_get(inv_Aup, jj, ii);
            sum += gsl_matrix_get(lapl_Adown, ii, jj) * gsl_matrix_get(inv_Adown, jj, ii);
        }
    }

    gsl_matrix_free(inv_Aup);
    gsl_matrix_free(inv_Adown);
    gsl_matrix_free(lapl_Aup);
    gsl_matrix_free(lapl_Adown);

    return -0.5 * sum;
}

double localEnergy2(gsl_matrix* Aup, gsl_matrix* Adown, double Rup[DIM][N2], double Rdown[DIM][N2], double r1[N2], double r2[N2], double eta){
     // Matrices of laplacians
    gsl_matrix* grad_Aup[DIM] = {gsl_matrix_alloc(N2, N2), gsl_matrix_alloc(N2, N2), gsl_matrix_alloc(N2, N2)};
    initialize_grad_mat(grad_Aup, Rup, r1, eta );
    gsl_matrix* grad_Adown[DIM] = {gsl_matrix_alloc(N2, N2), gsl_matrix_alloc(N2, N2), gsl_matrix_alloc(N2, N2)};
    initialize_grad_mat(grad_Adown, Rdown, r2, eta );

    double sum = 0.0, temp = 0.0;
    // Inverse matrices
    gsl_matrix* inv_Aup = invert_a_matrix(Aup);
    gsl_matrix* inv_Adown = invert_a_matrix(Adown);
    for(int kk = 0; kk < DIM; kk++)
    {
        for(int ii = 0; ii < N2; ii++)
        {
            for(int jj = 0; jj < N2; jj++)
            {
                for(int ll = 0; ll < N2; ll++)
                {
                    sum += gsl_matrix_get(grad_Aup[kk], ii, jj) * gsl_matrix_get(inv_Aup, jj, ii) * gsl_matrix_get(grad_Aup[kk], ii, ll) * gsl_matrix_get(inv_Aup, ll, ii);
                    sum += gsl_matrix_get(grad_Adown[kk], ii, jj) * gsl_matrix_get(inv_Adown, jj, ii)* gsl_matrix_get(grad_Adown[kk], ii, ll) * gsl_matrix_get(inv_Adown, ll, ii);
                } 
            }
        }
    }

    // for(int kk = 0; kk < DIM; kk ++) //i
    // {
    //     temp = 0;
    //     for(int ii = 0; ii < N2; ii++) // j
    //     {
    //         for(int jj = 0; jj < N2; jj++) // k
    //         {
    //             double a = gsl_matrix_get(grad_Aup[kk], ii, jj) * gsl_matrix_get(inv_Aup, jj, ii);
    //             double b = gsl_matrix_get(grad_Adown[kk], ii, jj)*gsl_matrix_get(inv_Adown, jj, ii);
    //             temp += a + b;
        
    //         }
    //         //sum += temp;
    //         //cout << sum << endl;
    //     }
    //     sum += temp*temp;

    // }

    gsl_matrix_free(inv_Aup);
    gsl_matrix_free(inv_Adown);
    for(int ii=0; ii < DIM; ii++)
    {
        gsl_matrix_free(grad_Aup[ii]);
        gsl_matrix_free(grad_Adown[ii]);
    }


    return 0.5 * sum;
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
// -----------------------------------------------------------------------------------


// int randomGenerator_int(int low, int up){
//     /* Generates random integer numbers in the range [low, up] */
//   return (rand() % (up - low + 1)) +low;
// }

// double randomGenerator_double(double low, double up){
//     double range = up - low;
//     double div = RAND_MAX / range;
//     return low + (rand()/div);
// }



// void gradient_of_matrix_x(gsl_matrix *matrix, double eta, double r[N2], double R[DIM][N2])
// {   /* Calculates one component of the gradient of the matrix */
//     double val = 0.0;

//         for (size_t ii = 0; ii < size; ii++) // run over rows : 4 positions of electrons
//         {
//                 val = - R[0][ii]/ (eta * eta) * Psi0(eta, r[ii]);
//                 gsl_matrix_set(matrix, ii, 0, val);
                
//                 val = ( 1.0 - R[0][ii] * R[0][ii] / (eta * eta) ) * Psi0(eta, r[ii]);
//                 gsl_matrix_set(matrix, ii, 1, val);

//                 val = - R[0][ii] * R[1][ii]/ (eta * eta) * Psi0(eta, r[ii]);
//                 gsl_matrix_set(matrix, ii, 2, val);
                            
//                 val = - R[0][ii] * R[2][ii]/ (eta * eta) * Psi0(eta, r[ii]);
//                 gsl_matrix_set(matrix, ii, 3, val);
                    
//         }
// }

// void gradient_of_matrix_y(gsl_matrix *matrix, double eta, double r[N2], double R[DIM][N2])
// {   /* Calculates one component of the gradient of the matrix */
//     double val = 0.0;

//         for (size_t ii = 0; ii < size; ii++) // run over rows : 4 positions of electrons
//         {
//                 val = - R[1][ii]/ (eta * eta) * Psi0(eta, r[ii]);
//                 gsl_matrix_set(matrix, ii, 0, val);
                
//                 val = - R[0][ii] * R[1][ii]/ (eta * eta) * Psi0(eta, r[ii]);
//                 gsl_matrix_set(matrix, ii, 1, val);

//                 val = ( 1.0 - R[1][ii] * R[1][ii] / (eta * eta) ) * Psi0(eta, r[ii]);
//                 gsl_matrix_set(matrix, ii, 2, val);
                            
//                 val = - R[1][ii] * R[2][ii]/ (eta * eta) * Psi0(eta, r[ii]);
//                 gsl_matrix_set(matrix, ii, 3, val);
                    
//         }
// }

// void gradient_of_matrix_z(gsl_matrix *matrix, double eta, double r[N2], double R[DIM][N2])
// {   /* Calculates one component of the gradient of the matrix */
//     double val = 0.0;

//         for (size_t ii = 0; ii < size; ii++) // run over rows : 4 positions of electrons
//         {
//                 val = - R[2][ii]/ (eta * eta) * Psi0(eta, r[ii]);
//                 gsl_matrix_set(matrix, ii, 0, val);
                
//                 val = - R[0][ii] * R[2][ii]/ (eta * eta) * Psi0(eta, r[ii]);
//                 gsl_matrix_set(matrix, ii, 1, val);

//                 val = - R[1][ii] * R[2][ii]/ (eta * eta) * Psi0(eta, r[ii]);
//                 gsl_matrix_set(matrix, ii, 2, val);
                            
//                 val = ( 1.0 - R[2][ii] * R[2][ii] / (eta * eta) ) * Psi0(eta, r[ii]);
//                 gsl_matrix_set(matrix, ii, 3, val);
                    
//         }
// }

// void laplacian_of_matrix(gsl_matrix *matrix, double eta, double r[N2], double R[DIM][N2])
// {   /* Calculated the laplacian of the matrix */
//     double val;
      
//         for (size_t ii = 0; ii < size; ii++)
//         {   val = 0.0;
//                 val = Psi0(eta, r[ii]) / (eta * eta) * (-3.0 + r[ii] * r[ii] / (eta * eta) );
//                 gsl_matrix_set(matrix, ii, 0, val);
//                 val = Psi0(eta, r[ii]) * R[0][ii] / (eta * eta) * (-5.0 + r[ii] * r[ii]/ (eta * eta) );
//                 gsl_matrix_set(matrix, ii, 1, val);
//                 val = Psi0(eta, r[ii]) * R[1][ii] / (eta * eta) * (-5.0 + r[ii] * r[ii]/ (eta * eta) );
//                 gsl_matrix_set(matrix, ii, 2, val);
//                 val = Psi0(eta, r[ii]) * R[2][ii] / (eta * eta) * (-5.0 + r[ii] * r[ii]/ (eta * eta) );
//                 gsl_matrix_set(matrix, ii, 3, val);
//         }
    
// }
    
// void GDtoDR(gsl_matrix *m1, gsl_matrix *m2, gsl_matrix *m3 , gsl_matrix *detInv, double GDratio[DIM][N2])
// {   /* Calculates the gradient determinant-to-determinant ratio when no move is performed */
//     double val = 0.0, inv = 0.0;
//     double totx= 0.0, toty = 0.0, totz = 0.0;
//     for (size_t ii = 0; ii < size; ii++) // run over rows : position of electrons
//     {  totx = toty = totz = 0.0;
//             for (size_t jj = 0; jj < size; jj++) // run over column : electrons
//             {  
//                 //  printf("val[%zu][%zu]: %lf \t inv[%zu][%zu]: %lf\t v*i: %lf\t tot: %lf\n", ii, jj, val, ii, jj, inv, val*inv, totz);
//                 val =  gsl_matrix_get(m1, ii, jj);  inv = gsl_matrix_get(detInv, jj, ii );
//                 totx += val *inv;
            
//                 val =  gsl_matrix_get(m2, ii, jj);
//                 toty += val *inv;
                
//                 val = gsl_matrix_get(m3, ii, jj);
//                 totz += val * inv;
                
//             }
//             GDratio[0][ii] = totx; // x-component
//             GDratio[1][ii] = toty; // y-component
//             GDratio[2][ii] = totz; // z-component
  
//     }

// }

// void LDtoDR(gsl_matrix *mm, gsl_matrix *detInv, double LDratio[N2])
// {   /* Calculated the laplacian determinant-to determinant ratio */
//     double val, inv, tot;
//     for (size_t ii = 0; ii < size; ii++) // run over rows
//     {   tot = 0.0;
//         for (size_t jj = 0; jj < size; jj++) // run over columns
//         {   val = gsl_matrix_get(mm, ii, jj);   inv = gsl_matrix_get(detInv, jj, ii);
//             tot += val*inv;
//         }
//     LDratio[ii] = tot;
//     }
    
// }



// double localEnergy2(double LDtoDRup[N2], double LDtoDRdown[N2], double GDtoDRup[DIM][N2], double GDtoDRdown[DIM][N2])
// { /* Second way to implement the calculation of the local kinetic energy */

//     double val1 = 0.0;
//     double val2[4];
//     double GD[3][4];
//     double sum = 0.0;
//     for (int ii = 0; ii < DIM; ii++)
//     {
//         val2[ii] = 0.0;
        
//     }
    
//     for (int ii = 0; ii < N2; ii++)
//     {
//              val1 += LDtoDRup[ii] +  LDtoDRdown[ii];
//     }

//     for (int jj = 0; jj < N2; jj++)
//     {  
//        for (int ii = 0; ii < DIM; ii++)
//        {
//             GD[ii][jj] = GDtoDRup[ii][jj] + GDtoDRdown[ii][jj];
//        }
//     }

//     for (int jj = 0; jj < N2; jj++)
//     {
//         for (int ii = 0; ii < DIM; ii++)
//         {
//             val2[jj] += GD[ii][jj] * GD[ii][jj];
//         }
        
//     }

//     for (int ii = 0; ii < N2; ii++)
//     {
//         sum += val2[ii];
//     }
    
    
 
//     return 1.0 / 4.0 * (- val1 + sum);
// }


// void getRatio(gsl_matrix *matrix, gsl_matrix *detInv, double ratio[N2] ){
   
//     double val = 0.0; 
//     double inv = 0.0;
//     double tot = 0.0;
//     for (size_t ii = 0; ii < size; ii++)
//     {   tot = 0.0;
//         for (size_t jj = 0; jj < size; jj++)
//         {
//             val = gsl_matrix_get(matrix, ii, jj);   inv = gsl_matrix_get(detInv, ii, jj);
//             tot += val*inv;
//             printf("val: %lf\t inv: %lf\t tot: %lf\n", val, inv, tot);
//         }
//         ratio[ii] = tot;
//     }
    
// }

