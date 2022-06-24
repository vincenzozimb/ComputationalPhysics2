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

/* ------------------------------- PARAMETERS OF THE SYSTEM ------------------------------- */
/* Specify the atom type. Choose "Na" for sodium or "K" for potassium */
char atom[] = "Na";  
double rs, Ray, rho;
    if(!strcmp(atom,"Na")){
        rs = 3.93;
    }else if(!strcmp(atom,"K")){
        rs = 4.86;
    }else{
        printf("\n\n\nERROR: Atom name wrong!\n\n\n");
        return 0;
    }
    Ray = rs * pow((double)N,1.0/3.0); // radius of the cluster (of its harmonic part)
    rho = 3.0 / (4.0 * M_PI * rs * rs * rs); // density of the jellium


/* Initial positions */
double Xold[DIM][N2]; // initial positions : 3D and 4 electrons
double Xnew[DIM][N2]; 
double L = 3.0 * Ray;
double h = L/(DIMx -1.0);
double x[DIMx], Vex[DIMx];

/* Set matrix of positions at 0 */
for (int ii = 0; ii < DIM; ii++)
{
    for (int jj = 0; jj < N2; jj++)
    {
        Xold[ii][jj] = Xnew[ii][jj] = 0.0;
    }
}

/* Create a mesh of positions and the external potential Vex */
for (int ii = 0; ii < DIMx; ii++)
{
    x[ii] = ii * h;
    Vex[ii] = Vext(x[ii],Ray, rho);
}

FILE *pf;
pf = fopen("dataVext.csv", "w");
fprint_two_vec(pf, x, Vex, DIMx);
fclose(pf);

double r2[4], r[4]; // radius squared and radius
for (int jj = 0; jj < 4; jj++)
{
    r2[jj] = Xold[0][jj] * Xold[0][jj] + Xold[1][jj] * Xold[1][jj] + Xold[2][jj] * Xold[2][jj]; 
    r[jj] = sqrt(r2[jj]);
    //printf("r2[%d] = %lf\n", jj, r2[jj]);
}


/* Variational parameters */
double eta[DIMe];
double eta_init = 0.2;
double eta_end = 10.5;
double deta = (eta_end - eta_init)/(DIMe);

/* Energy and variance */
double Etot_1 = 0.0, Etot_2 = 0.0;
double Etot2_1 = 0.0, Etot2_2 = 0.0;
double Var = 0.0;
double dE_1 = 0.0, dE_2 = 0.0;
double Energies_1[DIMe];
double Energies_2[DIMe];
double Variances[DIMe];

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

int index;

/* ------------------------------- EVOLUTION OF THE SYSTEM : 4 NON-INTERACTING ELECTRONS ------------------------------- */
for (int ee = 0; ee < DIMe; ee++)
{
    printf("\nciclo : %d #####################################################\n", ee);
    eta[ee] = eta_init + ee * deta; // initialize variational parameters

    Etot_1 = 0.0; // initialize total energy at each cycle
    Etot2_1 = 0.0;
    Etot_2 = 0.0; // initialize total energy at each cycle
    Etot2_2 = 0.0;

/* ########################################## Initialization ############################################# */
    /* Initialization of the matrix of positions : column -> electron, row -> x,y,z ( position in which
     I'm considering the electron) */
   // printf("Xold initial: \n");
   //printf("Xold\n");
    for (int ii = 0; ii < DIM; ii++)
    {
        for (int jj = 0; jj < N2; jj++)
        {
            Xold[ii][jj] = h * (rand()/(double)RAND_MAX - 1.0/2.0);
        }
    }

    /* Calculation of r */
    for (int ii = 0; ii < N2; ii++)
    {
        r2[ii] = Xold[0][ii] * Xold[0][ii] + Xold[1][ii] * Xold[1][ii] + Xold[2][ii] * Xold[2][ii]; 
        r[ii] = sqrt(r2[ii]);
    }
    
    /* Aup */
    // initialization_of_wf(Aup, gradAup_x, gradAup_y, gradAup_z, lapAup, GDratio_up, LDratio_up, detAup, eta[ee], r, Xold);
    initialize_mat_contents(Aup, eta[ee], r, Xold);
    // printf("Aup\n");
    // print_mat_contents(Aup);
    // printf("\n");
    gradient_of_matrix(gradAup_x, eta[ee], r, Xold, 0); 
    // printf("gradAup_x\n");
    // print_mat_contents(gradAup_x);
    // printf("\n");
    gradient_of_matrix(gradAup_y, eta[ee], r, Xold, 1); 
    // printf("gradAup_y\n");
    // print_mat_contents(gradAup_y);
    // printf("\n");
    gradient_of_matrix(gradAup_z, eta[ee], r, Xold, 2);
    // printf("gradAup_z\n");
    // print_mat_contents(gradAup_z);
    // printf("\n");
    laplacian_of_matrix(lapAup, eta[ee], r, Xold);
    // printf("lapAup\n");
    // print_mat_contents(lapAup);
    // printf("\n");
    // printf("\n");
    Aup_inv = invert_a_matrix(Aup); // To do because it doesn't return the inverse matrix
    // printf("Aup_inv\n");
    // print_mat_contents(Aup_inv);
    // printf("\n");
    initialize_mat_contents(Aup, eta[ee], r, Xold); // To do because invert_a_matrix reverts the matrix Aup
    GDtoDR_old(gradAup_x, gradAup_y, gradAup_z, Aup_inv, GDratio_up);
    // printf("GDtoDR\n"); 
    // for (int ii = 0; ii < DIM; ii++)
    // {
    //    for (int jj = 0; jj < N2; jj++)
    //    {
    //         printf("%lf\t", GDratio_up[ii][jj]);
    //    }
    //    printf("\n");
    // }
    // printf("\n");
    LDtoDR(lapAup, Aup_inv, LDratio_up); // LDtoD ratio
    // printf("LDtoDR\n"); 
    // for (int ii = 0; ii < N2; ii++)
    // {
    //     printf("%lf\n",LDratio_up[ii] );
    // }
    // printf("\n");
    detAup = determinant_of_matrix(Aup); // Slater determinant
    // printf("detAup : %lf\n\n", detAup);
    initialize_mat_contents(Aup, eta[ee], r, Xold); // Re-initialize m
    
    /* Adown */
    // initialization_of_wf(Adown, gradAdown_x, gradAdown_y, gradAdown_z, lapAdown, GDratio_down, LDratio_down, detAdown, eta[ee], r, Xold);
    initialize_mat_contents(Adown, eta[ee], r, Xold); //Initialize matrix Aup
    gradient_of_matrix(gradAdown_x, eta[ee], r, Xold, 0); // Gradient of Adown - Adown_x
    gradient_of_matrix(gradAdown_y, eta[ee], r, Xold, 1); // Gradient of Adown - Adown_y
    gradient_of_matrix(gradAdown_z, eta[ee], r, Xold, 2); // Gradient of Adown - Adown_z
    laplacian_of_matrix(lapAdown, eta[ee], r, Xold); // Laplacian of Adown
    Adown_inv = invert_a_matrix(Adown); // Inverse of Adown
    initialize_mat_contents(Adown, eta[ee], r, Xold); // Re-initialize Adown
    GDtoDR_old(gradAdown_x,gradAdown_y, gradAdown_z, Adown_inv, GDratio_down); //GDtoD ratio
    LDtoDR(lapAdown, Adown_inv, LDratio_down); // LDtoD ratio
    detAdown = determinant_of_matrix(Adown); // Slater determinant
    initialize_mat_contents(Adown, eta[ee], r, Xold); // Re-initialize Adown
    /* ################################################ Evolution ################################################# */

    psiOld = detAup * detAdown; // initial trial wave function with electrons in specific positions

    srand(time(0)); 
    
    /* ------------------------------------------------ Monte Carlo loop ------------------------------------------------*/
    for (int mm = 0; mm < M; mm++)
    {
        index = randomGenerator(0 , 3); // take the evolving electron in a random way
        //printf("index: %d\n", index);

        /* Update positions */
        for (int ii = 0; ii < DIM; ii++)
        {
            for (int jj = 0; jj < N2; jj++)
            {
                if (jj == index)
                {
                    Xnew[ii][jj] = Xold[ii][jj] + h * (rand()/(double)RAND_MAX - 1.0/2.0);
                } else {
                    Xnew[ii][jj] = Xold[ii][jj];
                }
                    
                }
    
            }    

            // for(int ii = 0; ii< DIM; ii++){
            //     for (int jj = 0; jj < N2; jj++)
            //     {
            //         printf("%lf\t", Xnew[ii][jj]);
            //     }
            //     printf("\n");
            // }
            // printf("\n");
            
        /*  Evaluation of wf in Xnew */
         // evolution_of_wf(Aup, Aup_inv, gradAup_x, gradAup_y, gradAup_z, lapAup, GDratio_up, LDratio_up, detAup, eta[ee], r, Xnew);
        initialize_mat_contents(Aup, eta[ee], r, Xnew); //Initialize matrix Aup
        gradient_of_matrix(gradAup_x, eta[ee], r, Xnew, 0); 
        gradient_of_matrix(gradAup_y, eta[ee], r, Xnew, 1); 
        gradient_of_matrix(gradAup_z, eta[ee], r, Xnew, 2); 
        laplacian_of_matrix(lapAup, eta[ee], r, Xnew); 
        Aup_inv = invert_a_matrix(Aup); 
        initialize_mat_contents(Aup, eta[ee], r, Xnew); 
        GDtoDR_new(Aup, gradAup_x,gradAup_y, gradAup_z, Aup_inv, GDratio_up); 
             printf("GDtoDR\n"); 
            for (int ii = 0; ii < DIM; ii++)
            {
                for (int jj = 0; jj < N2; jj++)
                {
                        printf("%lf\t", GDratio_up[ii][jj]);
                }
                printf("\n");
            }
             printf("\n");
        LDtoDR(lapAup, Aup_inv, LDratio_up); 
            printf("LDtoDR\n"); 
            for (int ii = 0; ii < N2; ii++)
            {
                printf("%lf\n",LDratio_up[ii] );
            }
            printf("\n");
        detAup = determinant_of_matrix(Aup); 
        initialize_mat_contents(Aup, eta[ee], r, Xnew); 
        initialize_mat_contents(Adown, eta[ee], r, Xnew); //Initialize matrix Aup
        gradient_of_matrix(gradAdown_x, eta[ee], r, Xnew, 0); // Gradient of Adown - Adown_x
        gradient_of_matrix(gradAdown_y, eta[ee], r, Xnew, 1); // Gradient of Adown - Adown_y
        gradient_of_matrix(gradAdown_z, eta[ee], r, Xnew, 2); // Gradient of Adown - Adown_z
        laplacian_of_matrix(lapAdown, eta[ee], r, Xnew); // Laplacian of Adown
        Adown_inv = invert_a_matrix(Adown); // Inverse of Adown
        initialize_mat_contents(Adown, eta[ee], r, Xnew); // Re-initialize Adown
        GDtoDR_new(Adown, gradAdown_x,gradAdown_y, gradAdown_z, Adown_inv, GDratio_down); 
        LDtoDR(lapAdown, Adown_inv, LDratio_down); // LDtoD ratio
        detAdown = determinant_of_matrix(Adown); // Slater determinant
        initialize_mat_contents(Adown, eta[ee], r, Xnew); // Re-initialize Adown
        
        psiNew = detAup * detAdown;
        // printf("psiNew: %lf\n", psiNew);
        /* Metropolis test */
        // printf("rand\t\t psiNew/psiOld\n");
        if (rand()/(double)RAND_MAX <= psiNew * psiNew/(psiOld * psiOld))
        {   
            // printf("%lf\t %lf\n", rand()/(double)RAND_MAX, psiNew * psiNew/(psiOld * psiOld) );
           for (int ii = 0; ii < DIM; ii++)
            {
                for (int jj = 0; jj < N2; jj++) // Update positions
                {
                    if (jj == index)
                    {
                        Xold[ii][jj] = Xnew[ii][jj];
                    } else {
                        Xold[ii][jj] = Xold[ii][jj];
                    }
                    
                }
    
            }    
            psiOld = psiNew; // Update wave function
            // printf("psiOld: %lf\n", psiOld);
            // printf("mm: %d\t 1\n", mm);
        // } else { printf("mm: %d\t 0\n", mm);
        }


        /* Compute local energy */
        // avoided thermalisation

        // control local energy anf GD and LD ratios (there should be a dependence on the old position )

        dE_1 = localEnergy1(LDratio_up, LDratio_down); 
        dE_2 = localEnergy2(LDratio_up, LDratio_down, GDratio_up, GDratio_down );
        Etot_1 += dE_1;
        Etot_2 += dE_2;
        Etot2_1 += dE_1 * dE_1;
        Etot2_2 += dE_2 * dE_2;
    
    }

    Etot_1 /= M;
    Etot2_1 /= M;
    Etot_2 /= M;
    Etot2_2 /= M;
    // Var = Etot2 - Etot *Etot;
    Energies_1[ee] = Etot_1;
    Energies_2[ee] = Etot_2;
    // Variances[ee] = Var;
}

/* Print results on file */
FILE *pf2;
pf2 = fopen("dataEnEl.csv", "w");
fprint_three_vec(pf2, eta, Energies_1,Energies_2, DIMe);
fclose(pf2);


/* Free memory space */
gsl_matrix_free(Aup);   gsl_matrix_free(Adown);
gsl_matrix_free(Aup_inv);   gsl_matrix_free(Adown_inv);
gsl_matrix_free(gradAup_x);   gsl_matrix_free(gradAup_y);   gsl_matrix_free(gradAup_z);
gsl_matrix_free(gradAdown_x);   gsl_matrix_free(gradAdown_y);   gsl_matrix_free(gradAdown_z);
gsl_matrix_free(lapAup);    gsl_matrix_free(lapAdown);

}