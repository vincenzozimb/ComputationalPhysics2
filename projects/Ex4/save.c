#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

#include "print_routines.h"
#define DIMx 1000
#define DIMe 100
#define N 8

static size_t size = 4;

/* ======================= FUNCTION HEADERS ======================= */
double Psi0(double eta, double r);
double Psi1(double eta, double r, double x);
void initialize_mat_contents(gsl_matrix *matrix, double eta, double r[4], double R[3][4]);
gsl_matrix * invert_a_matrix(gsl_matrix *matrix);
void print_mat_contents(gsl_matrix *matrix);
double determinant_of_matrix(gsl_matrix *matrix);


/* ======================= MAIN ======================= */
int main(){

/* Initial positions */
// double L = 8.0;
// double h = L / (DIMx - 1);
double X[3][4]; // initial positions : 3D and 4 electrons
double r2[4], r[4];
/* Initialize positions of 4 electrons - a caso */
X[0][0] = 1.0; X[1][0] = 0.2; X[2][0] = 0.5;
X[0][1] = 0.0; X[1][1] = 1.0; X[2][1] = 0.0; 
X[0][2]= 1.0; X[1][2] = 2.0; X[2][2] = 1.0;
X[0][3] = 0.5; X[1][3] = 1.0; X[2][3] = 1.5;

double posOld[3][4], posNew[3][4]; 
 
// for (int ii = 0; ii < 3; ii++)
// {
//     for (int jj = 0; jj < 4; jj++)
//     {
//         X[ii][jj] = rand()/(double)RAND_MAX;
//         printf("%lf\t", X[ii][jj]);
//     }
//    printf("\n");
// }
for (int jj = 0; jj < 4; jj++)
{
    r2[jj] = X[0][jj] * X[0][jj] + X[1][jj] * X[1][jj] + X[2][jj] * X[2][jj]; 
    r[jj] = sqrt(r2[jj]);
    //printf("r2[%d] = %lf\n", jj, r2[jj]);
}


/* Variational parameters */
double eta[DIMe];
double eta_init = 0.3;
double eta_end = 5.0;
double deta = (eta_end - eta_init)/(DIMe - 1.0);

/* Energy and variance */
double Etot = 0.0;
double Etot2 = 0.0;

/* Slater matrix */
gsl_matrix *Aup = gsl_matrix_alloc(size,size);

// Make 4 electrons evolve together
for (int ee = 0; ee < DIMe; ee++)
{
    eta[ii] = eta_init + ii * deta;

    Etot = 0.0;
    Etot2 = 0.0;


}
   initialize_mat_contents(Aup, eta[ee], r, X);
printf("Aup : \n");
print_mat_contents(Aup);
printf("\n");

// gsl_matrix_free(Aup);

/* Inverse of Aup matrix */
gsl_matrix *Aup_inv = invert_a_matrix(Aup);
printf("Inverse Aup: \n");
print_mat_contents(Aup_inv);
printf("\n");

/* Re-initialization*/
initialize_mat_contents(Aup, eta[45], r, X);
printf("Aup : \n");
print_mat_contents(Aup);
printf("\n");
printf("eta[]:%lf\n", eta[45]);

/* Slater determinant of A matrix */

double detAup;

detAup = determinant_of_matrix(Aup);
printf("\n\n detAup = %lf\n\n", detAup);
gsl_matrix_free(Aup);
}
/* ======================= FUNCTION BODIES ======================= */

double Psi0(double eta, double r){
    return exp(-r * r /(2.0 * eta * eta));
}

double Psi1(double eta, double r, double x){
    return x * exp(-r * r /(2.0 * eta * eta));
}

void initialize_mat_contents(gsl_matrix *matrix, double eta, double r[4], double R[3][4]){

    double val; 
    for ( size_t ii = 0; ii < size; ii++) // run over 4 positions of electrons
    {
        for ( size_t jj = 0; jj < size; jj++) // run over the wf of the 4 electrons
        {
            if(jj == 0){
                val = Psi0(eta, r[ii]);
                gsl_matrix_set(matrix, ii, jj, val);
            } else {
                 val = Psi1(eta, r[ii], R[ii][jj]);
                 gsl_matrix_set(matrix, ii, jj, val);
            }

        }
        
    }
    
}

gsl_matrix * invert_a_matrix(gsl_matrix *matrix)
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

void print_mat_contents(gsl_matrix *matrix)
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

double determinant_of_matrix(gsl_matrix *matrix){
int detMat;
    
    gsl_permutation *p = gsl_permutation_alloc(size);
    int s;
    int signum = 1;

    // Compute the LU decomposition of this matrix
    gsl_linalg_LU_decomp(matrix, p, &s);


    detMat = gsl_linalg_LU_det(matrix, signum);
    gsl_permutation_free(p);
return detMat;
}




    initialize_mat_contents(Adown, eta[ee], r, Xold); //Initialize matrix Aup
    gradient_of_matrix(gradAdown_x, eta[ee], r, Xold, 0); // Gradient of Adown - Adown_x
    gradient_of_matrix(gradAdown_y, eta[ee], r, Xold, 1); // Gradient of Adown - Adown_y
    gradient_of_matrix(gradAdown_z, eta[ee], r, Xold, 2); // Gradient of Adown - Adown_z
    laplacian_of_matrix(lapAdown, eta[ee], r, Xold); // Laplacian of Adown
    Adown_inv = invert_a_matrix(Adown); // Inverse of Adown
    initialize_mat_contents(Adown, eta[ee], r, Xold); // Re-initialize Adown
    GDtoDR_old(gradAdown_x,gradAdown_y, gradAdown_z, Adown_inv, GDratio_down); GDtoD ratio
    LDtoDR(lapAdown, Adown_inv, LDratio_down); // LDtoD ratio
    detAdown = determinant_of_matrix(Adown); // Slater determinant


    initialize_mat_contents(Aup, eta[ee], r, Xold); //Initialize matrix Aup
    // printf("Aup:\n");
    // print_mat_contents(Aup);
    // printf("\n");

    gradient_of_matrix(gradAup_x, eta[ee], r, Xold, 0); // Gradient of Aup - grAup_x
    // printf("gradAupx:\n");
    // print_mat_contents(gradAup_x);
    // printf("\n"); 
    gradient_of_matrix(gradAup_y, eta[ee], r, Xold, 1); // Gradient of Aup - grAup_y
    gradient_of_matrix(gradAup_z, eta[ee], r, Xold, 2); // Gradient of Aup - grAup_z

    laplacian_of_matrix(lapAup, eta[ee], r, Xold); // Laplacian of Aup
    // printf("lapAup: \n");
    // print_mat_contents(lapAup);
    // printf("\n");

    Aup_inv = invert_a_matrix(Aup); // Inverse of Aup

    initialize_mat_contents(Aup, eta[ee], r, Xold); // Re-initialize Aup

    GDtoDR_old(gradAup_x,gradAup_y, gradAup_z, Aup_inv, GDratio_up); // GDtoD ratio
    // printf("GDtoDR_old: \n"); 
    // print_mat_contents2(GDratio_up);
    // printf("\n");

    LDtoDR(lapAup, Aup_inv, LDratio_up); // LDtoD ratio
    // for (int ii = 0; ii < N2; ii++)
    // {
    //     printf("LDratio[%d] =%lf\t", ii, LDratio_up[ii]);
    // }
    // printf("\n");

    detAup = determinant_of_matrix(Aup); // Slater determinant
    //printf("det initial : %lf\n", detAup);



    printf("Adown\n");
    print_mat_contents(Adown);
    printf("\n");
    printf("gradAdown_x\n");
    print_mat_contents(gradAdown_x);
    printf("\n");
    printf("gradAup_y\n");
    print_mat_contents(gradAup_y);
    printf("\n");
    printf("gradAup_z\n");
    print_mat_contents(gradAup_z);
    printf("\n");
    printf("lapAdown\n");
    print_mat_contents(lapAdown);
    printf("\n");
    printf("\n");
    Adown_inv = invert_a_matrix(Adown); // To do because it doesn't return the inverse matrix
    printf("Adown_inv\n");
    print_mat_contents(Adown_inv);
    printf("\n");
    initialize_mat_contents(Adown, eta[ee], r, Xold); // To do because invert_a_matrix changes the matrix Aup
    printf("GDtoDR\n"); 
    for (int ii = 0; ii < DIM; ii++)
    {
       for (int jj = 0; jj < N2; jj++)
       {
            printf("%lf\t", GDratio_down[ii][jj]);
       }
       printf("\n");
    }
    printf("\n");
    printf("LDtoDR\n"); 
    for (int ii = 0; ii < N2; ii++)
    {
        printf("%lf\n",LDratio_down[ii] );
    }
    printf("\n");
    printf("detAdown : %lf\n\n", detAdown);
    initialize_mat_contents(Adown, eta[ee], r, Xold); // Re-initialize m