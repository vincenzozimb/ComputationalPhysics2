#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <string.h>
#include <time.h>

#include "util.h"
#include "print_routines.h"

int main(){
/* Initial and updated positions */
double Xold[DIM][N]; // 3D, 8 electrons
double Xold_up[DIM][N2], Xold_down[DIM][N2];
double Xnew_up[DIM][N2], Xnew_down[DIM][N2];
double Xnew[DIM][N];
double r_up[N2], r_down[N2];

/* Variational parameter */
double eta[DIMe];
double eta_init = 1.5;
double eta_end = 5.0;
double delta = 2.2 ; // parameter in the update of positions 
double deta = (eta_end - eta_init)/(DIMe);

/* Energy */
double Etot_1 = 0.0, Etot_2 = 0.0;
double Etot1_b = 0.0;
double Etot2_1 = 0.0, Etot2_2 = 0.0;
double dE_1 = 0.0, dE_2 = 0.0;
double dE1_b = 0.0;
double Energies_1[DIMe];
double Energies_1b[DIMe];
double Energies_2[DIMe];

/* Variance */
double Var_1 = 0.0, Var_2 = 0.0;
double Variances_1[DIMe];
double Variances_2[DIMe];

/* Statistical error */
double Error_1[DIMe];
double Error_2[DIMe];

/* Matrices */
gsl_matrix *Aup = gsl_matrix_alloc(size,size); // Slater matrix for 4 electrons with spin up
gsl_matrix *Aup_inv = gsl_matrix_alloc(size,size); // Inverse of the Slater matrix for 4 electrons with spin up
gsl_matrix *gradAup_x = gsl_matrix_alloc(size, size); // Component x of gradient of Aup
gsl_matrix *gradAup_y = gsl_matrix_alloc(size, size); // Component y of gradient of Aup
gsl_matrix *gradAup_z = gsl_matrix_alloc(size, size); // Component z of gradient of Aup
gsl_matrix *lapAup = gsl_matrix_alloc(size, size); // Laplacian of the matrix Aup
double detAup = 0.0; // Slater determinant of the matrix Aup
double GDratio_up[DIM][N2]; // Gradient determinant-to-determinant ratio for Aup
double LDratio_up[N2]; // Laplacian determinant-to-determinant ratio fo Aup

gsl_matrix *Adown = gsl_matrix_alloc(size,size); // Slater matrix for 4 electrons with spin down
gsl_matrix *Adown_inv = gsl_matrix_alloc(size,size); // Inverse of the Slater matrix for 4 electrons with spin down
gsl_matrix *gradAdown_x = gsl_matrix_alloc(size, size); // Component x of gradient of Adown
gsl_matrix *gradAdown_y = gsl_matrix_alloc(size, size); // Component y of gradient of Adown
gsl_matrix *gradAdown_z = gsl_matrix_alloc(size, size); // Component z of gradient of Adown
gsl_matrix *lapAdown = gsl_matrix_alloc(size, size); // Laplacian of the matrix Adown
double detAdown = 0.0; // Slater determinant of the matrix Adown
double GDratio_down[DIM][N2]; // Gradient determinant-to-determinant ratio for Adown
double LDratio_down[N2]; // Laplacian determinant-to-determinant ratio fo Adown

/* Total wave functions */
double psiOld, psiNew;
double ratio_up[N2]; // Rsd
double ratio_down[N2];

/* Indeces */
int cont = 0;

/* ratio over psi and random number */
double w, random;

srand(time(0)); // just call once to randomize the seed
/* ---------------------------------------------------------------------------------------------------------------------- */

/* Initial trial positions */
init_positions(Xold);
create_submat_and_r(Xold, Xold_up, Xold_down, r_up, r_down);

/* Loop over eta */
for (int ee = 0; ee < DIMe; ee++)
{  
    cont = 0;
    eta[ee] = eta_init + ee * deta;
    Etot_1 = 0.0; // initialize total energy at each cycle
    Etot2_1 = 0.0;
    Etot_2 = 0.0; // initialize total energy at each cycle
    Etot2_2 = 0.0;

    /* Initialize Slater matrices - Aup and Adown */
    init_Slater_matrices(Aup, Adown, eta[ee], r_up, r_down, Xold_up, Xold_down);

    /* Calculate Slater determinants */
    detAup = determinant_of_matrix(Aup); // Slater determinant
    detAdown = determinant_of_matrix(Adown); // Slater determinant
    init_Slater_matrices(Aup, Adown, eta[ee], r_up, r_down, Xold_up, Xold_down);

    /* Initial trial wave function */
    psiOld = detAup * detAdown; 

 /* ------------------------------------------------ Monte Carlo loop ------------------------------------------------*/
 for (int mm = 0; mm < M + thermal; mm++){
    /* Update the position of the moved electron */
    for (int index = 0; index < N; index++){// Loop over the 8 electrons
        for (int jj = 0; jj < N; jj++){
            for (int ii = 0; ii < DIM; ii++){
                if (jj == index){
                    Xnew[ii][jj] = Xold[ii][jj] + delta * (rand()/(double)RAND_MAX - 1.0/2.0);
                } else {
                    Xnew[ii][jj] = Xold[ii][jj];
                }
            }
        }
    create_submat_and_r(Xnew, Xnew_up, Xnew_down, r_up, r_down);
    /* Update wave function */
    init_Slater_matrices(Aup, Adown, eta[ee], r_up, r_down, Xnew_up, Xnew_down);
    detAup = determinant_of_matrix(Aup); // Slater determinant
    detAdown = determinant_of_matrix(Adown); // Slater determinant
    init_Slater_matrices(Aup, Adown, eta[ee], r_up, r_down, Xnew_up, Xnew_down); // re-initialize


    psiNew = detAup * detAdown;
    
    w = psiNew * psiNew/(psiOld * psiOld);
    random =  rand()/(double) RAND_MAX;
    if ( random <= w ){
        /* Accept the move of the electron */
        for (int jj = 0; jj < N; jj++)
        {
            for (int ii = 0; ii < DIM; ii++)
            {
               Xold[ii][jj] = Xnew[ii][jj];
            }
        }
    create_submat_and_r(Xold, Xold_up, Xold_down, r_up, r_down); 

    psiOld = psiNew;
    cont++;
    } // Closes Metropolis

    }// Closes loop over e

    /* Compute the local energy */
    if (mm >= thermal)
    {   
        /* Gradient of matrices Aup and Adown */
            /* Component x */
        gradient_of_matrix_x(gradAup_x, eta[ee], r_up, Xold_up); 
        gradient_of_matrix_x(gradAdown_x, eta[ee], r_down, Xold_down);
            /* Component y */
        gradient_of_matrix_y(gradAup_y, eta[ee], r_up, Xold_up);
        gradient_of_matrix_y(gradAdown_y, eta[ee], r_down, Xold_down);
            /* Component z */
        gradient_of_matrix_z(gradAup_z, eta[ee], r_up, Xold_up);
        gradient_of_matrix_z(gradAdown_z, eta[ee], r_down, Xold_down);
            /* Laplacian of matrices Aup and Adown */
        laplacian_of_matrix(lapAup, eta[ee], r_up, Xold_up);
        laplacian_of_matrix(lapAdown, eta[ee], r_down, Xold_down);
            /* Inverse of the matrices Aup and Adown */
        Aup_inv = invert_a_matrix(Aup);
        Adown_inv = invert_a_matrix(Adown);
            /* Re-initialize Aup and Adown - necessary because invert_a_matrix modifies the input matrix */
        init_Slater_matrices(Aup, Adown, eta[ee], r_up, r_down, Xold_up, Xold_down);

            /* Calculate ratios */
        GDtoDR( gradAup_x, gradAup_y, gradAup_z, Aup_inv, GDratio_up);
        GDtoDR( gradAdown_x, gradAdown_y, gradAdown_z, Adown_inv, GDratio_down);
      
        LDtoDR(lapAup, Aup_inv, LDratio_up); 
        LDtoDR(lapAdown, Adown_inv, LDratio_down); 

        dE1_b = localEnergy1_b(lapAup, lapAdown, Aup_inv, Adown_inv, gradAup_x, gradAup_y, gradAup_z, gradAdown_x, gradAdown_y, gradAdown_z);
        dE_1 = localEnergy1(LDratio_up, LDratio_down);
        dE_2 = localEnergy2(LDratio_up, LDratio_down, GDratio_up, GDratio_down );
        Etot_1 += dE_1;
        Etot1_b += dE1_b;
        Etot_2 += dE_2;
        Etot2_1 += dE_1 * dE_1;
        Etot2_2 += dE_2 * dE_2;
    }

} // Closes MC
printf("cont: %d\n", cont);
printf("Accept. prob: %lf\n", (double)cont/(8*(M+thermal)));
Etot_1 /= M;
Etot2_1 /= M;
Etot_2 /= M;
Etot2_2 /= M;
Etot1_b /= M;
Energies_1[ee] = Etot_1;
Energies_1b[ee] = Etot1_b;
Energies_2[ee] = Etot_2;


} // Closes loop over eta

/* Print results on file */
FILE *pf2, *pf3, *pf4;
pf2 = fopen("dataEnEl.csv", "w");
pf3 = fopen("dataVar.csv", "w");
pf4 = fopen("dataErr.csv", "w");
fprint_three_vec(pf2, eta, Energies_1,Energies_1b, DIMe);
fprint_three_vec(pf3, eta, Variances_1, Variances_2, DIMe);
fprint_three_vec(pf4, eta, Error_1, Error_2, DIMe);
fclose(pf2);
fclose(pf3);
fclose(pf4);



/* Free memory space */
gsl_matrix_free(Aup);   gsl_matrix_free(Adown);
gsl_matrix_free(Aup_inv);   gsl_matrix_free(Adown_inv);
gsl_matrix_free(gradAup_x);   gsl_matrix_free(gradAup_y);   gsl_matrix_free(gradAup_z);
gsl_matrix_free(gradAdown_x);   gsl_matrix_free(gradAdown_y);   gsl_matrix_free(gradAdown_z);
gsl_matrix_free(lapAup);    gsl_matrix_free(lapAdown);













}