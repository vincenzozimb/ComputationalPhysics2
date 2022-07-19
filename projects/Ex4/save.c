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
double eta_end = 4.0;
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

/* Define type of atom */
char atom[] = "Na"; // specify the atom type. Choose "Na" for sodium or "K" for potassium
double rs, R, rho;

    if(!strcmp(atom,"Na")){
        rs = 3.93;
    }else if(!strcmp(atom,"K")){
        rs = 4.86;
    }else{
        printf("\n\n\nERROR: Atom name wrong!\n\n\n");
        return 0;
    }
    R = rs * pow((double)N,1.0/3.0); // radius of the cluster (of its harmonic part)
    rho = 3.0 / (4.0 * M_PI * rs * rs * rs); // density of the jellium


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
//1 for (int index = 0; index < N; index++){// Loop over the 8 electrons
        for (int jj = 0; jj < N; jj++){
            for (int ii = 0; ii < DIM; ii++){
                //1 if (jj == index){
                    Xnew[ii][jj] = Xold[ii][jj] + delta * (rand()/(double)RAND_MAX - 1.0/2.0);
                //1} else {
                    //1Xnew[ii][jj] = Xold[ii][jj];
                //1}
            }
        }
    create_submat_and_r(Xnew, Xnew_up, Xnew_down, r_up, r_down);
    /* Update wave function */
    init_Slater_matrices(Aup, Adown, eta[ee], r_up, r_down, Xnew_up, Xnew_down);
    // printf("%d\tAup\n", index);
    // print_mat_contents(Aup);
    // printf("\n");
    // printf("%d\tAdown\n", index);
    // print_mat_contents(Adown);
    // printf("\n");

    detAup = determinant_of_matrix(Aup); // Slater determinant
    detAdown = determinant_of_matrix(Adown); // Slater determinant
    init_Slater_matrices(Aup, Adown, eta[ee], r_up, r_down, Xnew_up, Xnew_down); // re-initialize

    psiNew = detAup * detAdown;
    
    w = pow(psiNew/psiOld, 2.0);
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
    // printf("cont: %d\n", cont);
    } // Closes Metropolis

    //1}// Closes loop over e

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

        dE1_b = localEnergy1_b(LDratio_up, LDratio_down, Aup_inv, Adown_inv, gradAup_x, gradAup_y, gradAup_z, gradAdown_x, gradAdown_y, gradAdown_z);
        dE_1 = localEnergy1(LDratio_up, LDratio_down);
        dE_2 = localEnergy2(LDratio_up, LDratio_down, GDratio_up, GDratio_down );
        Etot_1 += dE_1 ;//+ Vext(r_up, r_down, R, rho);
        Etot1_b += dE1_b; //+ Vext(r_up, r_down, R, rho);
        Etot_2 += dE_2; //+ Vext(r_up, r_down, R, rho);
        // printf("dE1: %lf\t dE1_b: %lf\t dE2: %lf\t Etot1: %lf\t Etot1_b: %lf\t Etot2: %lf\n", dE_1, dE_2)
        Etot2_1 += dE_1 * dE_1;
        Etot2_2 += dE_2 * dE_2;
    }
//printf("%d\t %lf\t %lf\t%lf\n", mm, dE_1, dE_2, dE1_b);

} // Closes MC
printf("cont: %d\n", cont);
printf("Accept. prob: %lf\n", (double)cont/(8*(M+thermal)));
Etot_1 /= M;
Etot2_1 /= M;
Etot_2 /= M;
Etot2_2 /= M;
Etot1_b /= M;
Energies_1[ee] = Etot_1 ;
Energies_1b[ee] = Etot1_b;
Energies_2[ee] = Etot_2;
} // Closes loop over eta

/* Print results on file */
FILE *pf2, *pf3, *pf4;
pf2 = fopen("dataEnEl.csv", "w");
pf3 = fopen("dataVar.csv", "w");
pf4 = fopen("dataErr.csv", "w");
fprint_three_vec(pf2, eta, Energies_1b, Energies_1, DIMe);
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
    // printf("LDtoDR_down :\n"); 
    // for (int ii = 0; ii < N2; ii++)
    // {
    //     printf("%lf\n",LDratio_down[ii] );
    // }
    // printf("\n");

}

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
double eta_end = 4.0;
double delta = 1.2 ; // parameter in the update of positions 
double deta = (eta_end - eta_init)/(DIMe);

/* Energy */
double Etot_1 = 0.0, Etot_2 = 0.0;
double Etot2_1 = 0.0, Etot2_2 = 0.0;
double dE_1 = 0.0, dE_2 = 0.0;
double Energies_1[DIMe];
double Energies_2[DIMe];



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

/* Define type of atom */
char atom[] = "Na"; // specify the atom type. Choose "Na" for sodium or "K" for potassium
double rs, R, rho;

    if(!strcmp(atom,"Na")){
        rs = 3.93;
    }else if(!strcmp(atom,"K")){
        rs = 4.86;
    }else{
        printf("\n\n\nERROR: Atom name wrong!\n\n\n");
        return 0;
    }
    R = rs * pow((double)N,1.0/3.0); // radius of the cluster (of its harmonic part)
    rho = 3.0 / (4.0 * M_PI * rs * rs * rs); // density of the jellium


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
    // Move all electrons
        for (int jj = 0; jj < N; jj++){
            for (int ii = 0; ii < DIM; ii++){
                    Xnew[ii][jj] = Xold[ii][jj] + delta * (rand()/(double)RAND_MAX - 1.0/2.0);
            }
        }
    create_submat_and_r(Xnew, Xnew_up, Xnew_down, r_up, r_down);
    /* Update wave function */
    init_Slater_matrices(Aup, Adown, eta[ee], r_up, r_down, Xnew_up, Xnew_down);

    detAup = determinant_of_matrix(Aup); // Slater determinant
    detAdown = determinant_of_matrix(Adown); // Slater determinant
    init_Slater_matrices(Aup, Adown, eta[ee], r_up, r_down, Xnew_up, Xnew_down); // re-initialize

    psiNew = detAup * detAdown;
    
    w = pow(psiNew/psiOld, 2.0);
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

        dE_1 = localEnergy1(LDratio_up, LDratio_down);
        dE_2 = localEnergy2(LDratio_up, LDratio_down, GDratio_up, GDratio_down );
        Etot_1 += dE_1 + Vext(r_up, r_down, R, rho);
        Etot_2 += dE_2 + Vext(r_up, r_down, R, rho);
        Etot2_1 += dE_1 * dE_1;
        Etot2_2 += dE_2 * dE_2;
    }

} // Closes MC
printf("cont: %d\n", cont);
printf("Accept. prob: %lf\n", (double)cont/((M+thermal)));
Etot_1 /= M;
Etot2_1 /= M;
Etot_2 /= M;
Etot2_2 /= M;

Energies_1[ee] = Etot_1 ;
Energies_2[ee] = Etot_2;
} // Closes loop over eta

/* Print results on file */
FILE *pf2;
pf2 = fopen("dataEnEl.csv", "w");
fprint_three_vec(pf2, eta, Energies_1, Energies_2, DIMe);
fclose(pf2);


/* Free memory space */
gsl_matrix_free(Aup);   gsl_matrix_free(Adown);
gsl_matrix_free(Aup_inv);   gsl_matrix_free(Adown_inv);
gsl_matrix_free(gradAup_x);   gsl_matrix_free(gradAup_y);   gsl_matrix_free(gradAup_z);
gsl_matrix_free(gradAdown_x);   gsl_matrix_free(gradAdown_y);   gsl_matrix_free(gradAdown_z);
gsl_matrix_free(lapAup);    gsl_matrix_free(lapAdown);


}








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

void init_positions(double mat[DIM][N]){
    /* Initialize matrix by hand */
    // mat[0][0] = 0.5;   mat[1][0] = 1.0; mat[2][0] = 0.0;
    // mat[0][1] = 0.0;   mat[1][1] = 1.0; mat[2][1] = 1.0;
    // mat[0][2] = 1.0;   mat[1][2] = 1.0; mat[2][2] = 1.5;
    // mat[0][3] = 1.0;   mat[1][3] = 1.5; mat[2][3] = 1.0;
    // mat[0][4] = 0.0;   mat[1][4] = 0.5; mat[2][4] = 0.5;
    // mat[0][5] = 1.5;   mat[1][5] = 1.0; mat[2][5] = 0.0;
    // mat[0][6] = 0.0;   mat[1][6] = 1.0; mat[2][6] = 0.0;
    // mat[0][7] = 0.5;   mat[1][7] = 1.0; mat[2][7] = 1.5;

    for (int ii = 0; ii < DIM; ii++)
    {
        for (int jj = 0; jj < N; jj++)
        {
            mat[ii][jj] = (ii+1 + jj+1)/10.0;
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

double orbital(double X[DIM][N2], double r[N2], double eta, double ii, double jj){
// jj -> column : orbital / electron
// ii -> row : electron position
    if (jj == 0)
    {
        return Psi0(eta, r[ii]);
    } else {
        return PSi1(eta, r[ii], X[jj-1][ii]);
    }
    
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
            printf("%e ", element);
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
    for (size_t ii = 0; ii < size; ii++) // run over rows : position of electrons
    {  totx = toty = totz = 0.0;
            for (size_t jj = 0; jj < size; jj++) // run over column : electrons
            {  
                //  printf("val[%zu][%zu]: %lf \t inv[%zu][%zu]: %lf\t v*i: %lf\t tot: %lf\n", ii, jj, val, ii, jj, inv, val*inv, totz);
                val =  gsl_matrix_get(m1, ii, jj);  inv = gsl_matrix_get(detInv, jj, ii );
                totx += val *inv;
            
                val =  gsl_matrix_get(m2, ii, jj);
                toty += val *inv;
                
                val = gsl_matrix_get(m3, ii, jj);
                totz += val * inv;
                
            }
            GDratio[0][ii] = totx; // x-component
            GDratio[1][ii] = toty; // y-component
            GDratio[2][ii] = totz; // z-component
  
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


double localEnergy2(double LDtoDRup[N2], double LDtoDRdown[N2], double GDtoDRup[DIM][N2], double GDtoDRdown[DIM][N2])
{ /* Second way to implement the calculation of the local kinetic energy */

    double val1 = 0.0;
    double val2[4];
    double GD[3][4];
    double sum = 0.0;
    for (int ii = 0; ii < DIM; ii++)
    {
        val2[ii] = 0.0;
        
    }
    
    for (int ii = 0; ii < N2; ii++)
    {
             val1 += LDtoDRup[ii] +  LDtoDRdown[ii];
    }

    for (int jj = 0; jj < N2; jj++)
    {  
       for (int ii = 0; ii < DIM; ii++)
       {
            GD[ii][jj] = GDtoDRup[ii][jj] + GDtoDRdown[ii][jj];
       }
    }

    for (int jj = 0; jj < N2; jj++)
    {
        for (int ii = 0; ii < DIM; ii++)
        {
            val2[jj] += GD[ii][jj] * GD[ii][jj];
        }
        
    }

    for (int ii = 0; ii < N2; ii++)
    {
        sum += val2[ii];
    }
    
    
 
    return 1.0 / 4.0 * (- val1 + sum);
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

