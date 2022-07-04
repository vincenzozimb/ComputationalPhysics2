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
/* ------------------------------------------- PARAMETERS OF THE SYSTEM ------------------------------------------- */

/* Initial and updated positions */
double Xold[DIM][N]; // 3D, 8 electrons
double Xold_up[DIM][N2], Xold_down[DIM][N2];
double Xnew_up[DIM][N2], Xnew_down[DIM][N2];
double Xnew[DIM][N];

double r_up[N2], r_down[N2], r2_up[N2], r2_down[N2];

/* Variational parameters */
double eta[DIMe];
double eta_init = 2.0;
double eta_end = 6.0;
double deta = (eta_end - eta_init)/(DIMe);

/* Energy */
double Etot_1 = 0.0, Etot_2 = 0.0;
double Etot2_1 = 0.0, Etot2_2 = 0.0;
double dE_1 = 0.0, dE_2 = 0.0;
double Energies_1[DIMe];
double Energies_2[DIMe];

/* Variance */
double Var_1 = 0.0, Var_2 = 0.0;
double Variances_1[DIMe];
double Variances_2[DIMe];
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

/* Indeces */
int cont = 0;
int b = 0;
int index;

double delta = 1.5 ; // parameter in the update of positions 

srand(time(0)); // just call once to randomize the seed

/* Initial trial positions */
    for (int jj = 0; jj < N; jj++) // columns : run over electrons
    {
        for (int ii = 0; ii < DIM; ii++) // rows : run over components (x,y,z)
        {
            Xold[ii][jj] = randomGenerator_int(0, delta); // - delta * (rand()/(double)RAND_MAX - 1.0/2.0);
            if (jj < 4)
            {
                Xold_up[ii][jj] = Xold[ii][jj];
            } else {
                Xold_down[ii][jj-4] = Xold[ii][jj];
            }
            
        }
    }
    // printf("Xold: \n");
    // for (int ii = 0; ii < DIM; ii++)
    // {
    //     for (int jj = 0; jj < N; jj++)
    //     {
    //         printf("%lf\t", Xold[ii][jj]);
    //     }
    //     printf("\n");
    // }
    

    for (int jj = 0; jj < N2; jj++)
    {
        r2_up[jj] = Xold_up[0][jj] * Xold_up[0][jj] + Xold_up[1][jj] * Xold_up[1][jj] + Xold_up[2][jj] * Xold_up[2][jj];
        r_up[jj] = sqrt(r2_up[jj]);

        r2_down[jj] = Xold_down[0][jj] * Xold_down[0][jj] + Xold_down[1][jj] * Xold_down[1][jj] + Xold_down[2][jj] * Xold_down[2][jj];
        r_down[jj] = sqrt(r2_down[jj]);
        // printf("r_up[%d]: %lf\n", jj, r_up[jj]);
    }

/* ------------------------------- EVOLUTION OF THE SYSTEM : 8 NON-INTERACTING ELECTRONS ------------------------------- */
for (int ee = 0; ee < DIMe; ee++)
{   
    cont = 0;
    // printf("\nCYCLE : %d #########################################\n", ee);
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
    // printf("psiOld: %lf\t", psiOld);
  
    /* ------------------------------------------------ Monte Carlo loop ------------------------------------------------*/
for (int mm = 0; mm < M + thermal; mm++){
    b = 0;
    /* Select one of the 8 electrons randomly */
     index = randomGenerator_int(0 , N-1); // take the evolving electron in a random way 

    /* Update the position of that electron and r */
    
    for (int jj = 0; jj< N; jj++)
    {
        for (int ii = 0; ii < DIM; ii++)
        {
            if (jj == index)
            {
                Xnew[ii][jj] = Xold[ii][jj] + delta * (rand()/(double)RAND_MAX - 1.0/2.0);
            } else {
                Xnew[ii][jj] = Xold[ii][jj];
            }
            if (jj < 4)
            {
                Xnew_up[ii][jj] = Xnew[ii][jj];
            } else {
                Xnew_down[ii][jj-4] = Xnew[ii][jj];
            }
        
        }
        
    }
    // printf("Xnew: \n");
    // for (int ii = 0; ii < DIM; ii++)
    // {
    //     for (int jj = 0; jj < N; jj++)
    //     {
    //         printf("%lf\t", Xnew[ii][jj]);
    //     }
    //     printf("\n");
    // }
    for (int jj = 0; jj < N2; jj++)
    {
        r2_up[jj] = Xnew_up[0][jj] * Xnew_up[0][jj] + Xnew_up[1][jj] * Xnew_up[1][jj] + Xnew_up[2][jj] * Xnew_up[2][jj];
        r_up[jj] = sqrt(r2_up[jj]);

        r2_down[jj] = Xnew_down[0][jj] * Xnew_down[0][jj] + Xnew_down[1][jj] * Xnew_down[1][jj] + Xnew_down[2][jj] * Xnew_down[2][jj];
        r_down[jj] = sqrt(r2_down[jj]);
    }

    /* Update wave function */
    init_Slater_matrices(Aup, Adown, eta[ee], r_up, r_down, Xnew_up, Xnew_down);
    detAup = determinant_of_matrix(Aup); // Slater determinant
    detAdown = determinant_of_matrix(Adown); // Slater determinant
    init_Slater_matrices(Aup, Adown, eta[ee], r_up, r_down, Xnew_up, Xnew_down); // re-initialize

    psiNew = detAup * detAdown;
    //printf("% d\t psiNew: %lf\n", mm ,psiNew);
    
    // printf("dAu:%lf\t dAd: %lf\t dAuO: %lf\t dAdO: %lf\t Rsd: %lf\n", R_sd, detAup, detAdown, detAup_old, detAdown_old);
    /* Metropolis test*/
    double w, random; // acceptance ratio
    w = psiNew * psiNew/(psiOld * psiOld);
    random =  rand()/(double) RAND_MAX;
    if ( random <= w )
    {
         // printf("rand: %lf\t psiN/psiO: %lf\t", random, psiNew * psiNew/(psiOld * psiOld) );
        /* Update positions (new become old) and r */ 
       
        for (int jj = 0; jj < N; jj++)
        {
            for (int ii = 0; ii < DIM; ii++)
            {
               Xold[ii][jj] = Xnew[ii][jj];
               if (jj < 4)
                {
                    Xold_up[ii][jj] = Xold[ii][jj];
                } else {
                    Xold_down[ii][jj-4] = Xold[ii][jj];
                }
            }
            
        }
        for (int jj = 0; jj < N2; jj++)
        {
            r2_up[jj] = Xold_up[0][jj] * Xold_up[0][jj] + Xold_up[1][jj] * Xold_up[1][jj] + Xold_up[2][jj] * Xold_up[2][jj];
            r_up[jj] = sqrt(r2_up[jj]);

            r2_down[jj] = Xold_down[0][jj] * Xold_down[0][jj] + Xold_down[1][jj] * Xold_down[1][jj] + Xold_down[2][jj] * Xold_down[2][jj];
            r_down[jj] = sqrt(r2_down[jj]);
        }

        /* Update wave function */
        psiOld = psiNew;  
        cont++;       
        b = 1;
    }

    /* Compute the local energy */
    if (mm > thermal)
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

        dE_1 = localEnergy1(LDratio_up, LDratio_down);
        // dE_1 = -1.0/2.0 * (LDratio_up[0] + LDratio_down[0] + LDratio_up[1] + LDratio_down[1] + LDratio_up[2] + LDratio_down[2] + LDratio_up[3] + LDratio_down[3] );
        dE_2 = localEnergy2(LDratio_up, LDratio_down, GDratio_up, GDratio_down );
        // printf("%d\tdE2: %lf\n", mm ,dE_2);
        Etot_1 += dE_1;
        Etot_2 += dE_2;
        Etot2_1 += dE_1 * dE_1;
        Etot2_2 += dE_2 * dE_2;

    }
}
 //printf("cont: %d\n", cont-thermal);
 // printf("Accept. prob: %lf\n", (double)cont/(M+thermal));
Etot_1 /= M;
Etot2_1 /= M;
Etot_2 /= M;
Etot2_2 /= M;
// printf("En1: %lf \t En2: %lf\n", Etot_1,Etot_2);
Var_1 = Etot2_1 - Etot_1 * Etot_1;
Var_2 = Etot2_2 - Etot_2 * Etot_2;
Energies_1[ee] = Etot_1;
Energies_2[ee] = Etot_2;
Variances_1[ee] = Var_1;
Variances_2[ee] = Var_2;
Error_1[ee] = sqrt(1.0/(M-1.0) *(Etot2_1 - Etot_1 * Etot_1));
Error_2[ee] = sqrt(1.0/(M-1.0) *(Etot2_2 - Etot_2 * Etot_2));
}


/* Print results on file */
FILE *pf2, *pf3, *pf4;
pf2 = fopen("dataEnEl.csv", "w");
pf3 = fopen("dataVar.csv", "w");
pf4 = fopen("dataErr.csv", "w");
fprint_three_vec(pf2, eta, Energies_1,Energies_2, DIMe);
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

// /* Gradient of matrices Aup and Adown */
// /* Component x */
//     gradient_of_matrix_x(gradAup_x, eta[ee], r_up, Xold_up); 
//     printf("gradAup_x\n");
//     print_mat_contents(gradAup_x);
//     printf("\n");
//     gradient_of_matrix_x(gradAdown_x, eta[ee], r_down, Xold_down);
//     printf("gradAdown_x\n");
//     print_mat_contents(gradAdown_x);
//     printf("\n");
// /* Component y */
//     gradient_of_matrix_y(gradAup_y, eta[ee], r_up, Xold_up); 
//     printf("gradAup_y\n");
//     print_mat_contents(gradAup_y);
//     printf("\n");
//     gradient_of_matrix_y(gradAdown_y, eta[ee], r_down, Xold_down); 
//     printf("gradAdown_y\n");
//     print_mat_contents(gradAdown_y);
//     printf("\n");
// /* Component z */
//     gradient_of_matrix_z(gradAup_z, eta[ee], r_up, Xold);
//     printf("gradAup_z\n");
//     print_mat_contents(gradAup_z);
//     printf("\n");
//     gradient_of_matrix_z(gradAdown_z, eta[ee], r_down, Xold_down); 
//     printf("gradAdown_z\n");
//     print_mat_contents(gradAdown_z);
//     printf("\n");

// /* Laplacian of matrices Aup and Adown */
//     laplacian_of_matrix(lapAup, eta[ee], r_up, Xold_up);
//     printf("lapAup\n");
//     print_mat_contents(lapAup);
//     printf("\n");
//     laplacian_of_matrix(lapAdown, eta[ee], r_down, Xold_down);
//     printf("lapAdown\n");
//     print_mat_contents(lapAdown);
//     printf("\n");
// /* Inverse of the matrices Aup and Adown */
//     Aup_inv = invert_a_matrix(Aup);
//     printf("Aup_inv\n");
//     print_mat_contents(Aup_inv);
//     printf("\n");
//     Adown_inv = invert_a_matrix(Adown);
//     printf("Adown_inv\n");
//     print_mat_contents(Adown_inv);
//     printf("\n");
// /* Re-initialize Aup and Adown - necessary because invert_a_matrix modifies the input matrix */
//     initialize_mat_contents(Aup, eta[ee], r_up, Xold_up); 
//     initialize_mat_contents(Adown, eta[ee], r_down, Xold_down); 
// /* Calculate ratios */
//     GDtoDR_old(gradAup_x, gradAup_y, gradAup_z, Aup_inv, GDratio_up);
//     printf("GDtoDR_up: \n"); 
//     for (int ii = 0; ii < DIM; ii++)
//     {
//        for (int jj = 0; jj < N2; jj++)
//        {
//             printf("%lf\t", GDratio_up[ii][jj]);
//        }
//        printf("\n");
//     }
//     printf("\n");
//     GDtoDR_old(gradAdown_x, gradAdown_y, gradAdown_z, Adown_inv, GDratio_down);
//     printf("GDtoDR_down\n"); 
//     for (int ii = 0; ii < DIM; ii++)
//     {
//        for (int jj = 0; jj < N2; jj++)
//        {
//             printf("%lf\t", GDratio_down[ii][jj]);
//        }
//        printf("\n");
//     }
//     printf("\n");
//   LDtoDR(lapAup, Aup_inv, LDratio_up); // LDtoD ratio
//     printf("LDtoDR_up :\n"); 
//     for (int ii = 0; ii < N2; ii++)
//     {
//         printf("%lf\n",LDratio_up[ii] );
//     }
//     printf("\n");
//     LDtoDR(lapAdown, Adown_inv, LDratio_down); // LDtoD ratio
//     printf("LDtoDR_down :\n"); 
//     for (int ii = 0; ii < N2; ii++)
//     {
//         printf("%lf\n",LDratio_down[ii] );
//     }
//     printf("\n");

}


for (int jj = 0; jj< N; jj++)
    {
        for (int ii = 0; ii < DIM; ii++)
        {
        
            Xnew[ii][jj] = Xold[ii][jj] + delta * (rand()/(double)RAND_MAX - 1.0/2.0);
            if (jj < 4)
            {
                Xnew_up[ii][jj] = Xnew[ii][jj];
            } else {
                Xnew_down[ii][jj-4] = Xnew[ii][jj];
            }
        
        }
        
    }



    double localEnergy1( double LDtoDRup[N2], double LDtoDRdown[N2]){
    /* First way to implement the calculation of the local kinetic energy */
   double val = 0.0;
    for (int ii = 0; ii < N2; ii++)
    {
        val += LDtoDRup[ii] + LDtoDRdown[ii] ;
    }
    return -1.0/2.0 * val;
}