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
double eta_end = 1.5;
double deta = (eta_end - eta_init)/(DIMe);

/* Energy and variance */
double Etot_1 = 0.0;
double Etot_2 = 0.0;
double Etot2_1 = 0.0;
double Etot2_2 = 0.0;
double Var = 0.0;
double dE_1, dE_2;
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
            //printf("%lf\t", Xold[ii][jj]);
        }
        //printf("\n");
    }
    //printf("\n");

    /* Calculation of r */
    for (int ii = 0; ii < N2; ii++)
    {
        r2[ii] = Xold[0][ii] * Xold[0][ii] + Xold[1][ii] * Xold[1][ii] + Xold[2][ii] * Xold[2][ii]; 
        r[ii] = sqrt(r2[ii]);
    }
    
    /* Aup */
    initialization_of_wf(Aup, Aup_inv, gradAup_x, gradAup_y, gradAup_z, lapAup, GDratio_up, LDratio_up, detAup, eta[ee], r, Xold);
    // printf("Aup\n");
    // print_mat_contents(Aup);
    // printf("\n");
    // printf("GDratio_up\n");
    // for (int ii = 0; ii < DIM; ii++)
    // {
    //     for (int jj = 0; jj < N2; jj++)
    //     {
    //        printf("%lf\t", GDratio_up[ii][jj]);
    //     }
    //     printf("\n");
        
    // }
    
    
    /* Adown */
    initialization_of_wf(Adown, Adown_inv, gradAdown_x, gradAdown_y, gradAdown_z, lapAdown, GDratio_down, LDratio_down, detAdown, eta[ee], r, Xold);

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
            
        /*  Evaluation of wf in Xnew */
        evolution_of_wf(Aup, Aup_inv, gradAup_x, gradAup_y, gradAup_z, lapAup, GDratio_up, LDratio_up, detAup, eta[ee], r, Xnew);
        evolution_of_wf(Adown, Adown_inv, gradAdown_x, gradAdown_y, gradAdown_z, lapAdown, GDratio_down, LDratio_down, detAdown, eta[ee], r, Xold);

        psiNew = detAup * detAdown;

    
        /* Metropolis test */
        if (rand()/(double)RAND_MAX <= psiNew * psiNew/(psiOld * psiOld))
        {
           for (int ii = 0; ii < DIM; ii++)
            {
                for (int jj = 0; jj < N2; jj++)
                {
                    if (jj == index)
                    {
                        Xold[ii][jj] = Xnew[ii][jj];
                    } else {
                        Xold[ii][jj] = Xold[ii][jj];
                    }
                    
                }
    
            }    
            psiOld = psiNew; // Updated wave function
        }

        /* Compute local energy */
        // avoided thermalisation
        dE_1 = localEnergy1( LDratio_up, LDratio_down); 
        dE_2 = localEnergy2( LDratio_up, LDratio_down, GDratio_up, GDratio_down );
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