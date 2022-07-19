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
double eta_init = 1.6;
double eta_end = 4.0;
double delta = 1.4; // parameter in the update of positions 
double deta = (eta_end - eta_init)/(DIMe);

/* Energy */
double Etot_1 = 0.0, Etot_2 = 0.0;
double Etot2_1 = 0.0, Etot2_2 = 0.0;
double dE_1 = 0.0, dE_2 = 0.0;
double Energies_1[DIMe];
double Energies_2[DIMe];


/* Matrices */
// gsl_matrix *Aup = gsl_matrix_alloc(size,size); // Slater matrix for 4 electrons with spin up
gsl_matrix *Aup_inv = gsl_matrix_alloc(size,size); // Inverse of the Slater matrix for 4 electrons with spin up
gsl_matrix *gradAup_x = gsl_matrix_alloc(size, size); // Component x of gradient of Aup
gsl_matrix *gradAup_y = gsl_matrix_alloc(size, size); // Component y of gradient of Aup
gsl_matrix *gradAup_z = gsl_matrix_alloc(size, size); // Component z of gradient of Aup
gsl_matrix *lapAup = gsl_matrix_alloc(size, size); // Laplacian of the matrix Aup
double detAup = 0.0; // Slater determinant of the matrix Aup
double GDratio_up[DIM][N2]; // Gradient determinant-to-determinant ratio for Aup
double LDratio_up[N2]; // Laplacian determinant-to-determinant ratio fo Aup

// gsl_matrix *Adown = gsl_matrix_alloc(size,size); // Slater matrix for 4 electrons with spin down
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

double mmm[M]; double Energies_1m[M]; double dE_1m[M];

srand(time(0)); // just call once to randomize the seed (time0)

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

for (int ee = 0; ee < DIMe; ee++){
gsl_matrix *Aup = gsl_matrix_alloc(size,size);
gsl_matrix *Adown = gsl_matrix_alloc(size,size);
    Etot_1 = 0.0; Etot_2 = 0.0; 
    eta[ee] = eta_init + ee * deta;
    cont = 0;
    // initial position 
    init_positions(Xold, delta);
    create_submat_and_r(Xold, Xold_up, Xold_down, r_up, r_down);
    init_Slater_matrices(Aup, Adown, eta[ee], r_up, r_down, Xold_up, Xold_down);
    /* Calculate Slater determinants */
    detAup = determinant_of_matrix(Aup); // Slater determinant
    detAdown = determinant_of_matrix(Adown); // Slater determinant
    init_Slater_matrices(Aup, Adown, eta[ee], r_up, r_down, Xold_up, Xold_down);
    psiOld = detAup * detAdown; 

    for (int mm = 0; mm < M + thermal; mm++){ // MC
    if (mm >= thermal){ // calculation of energy
    dE_1 = localEnergy1(Aup, Adown, Xold_up, Xold_down, r_up, r_down, eta[ee]);
    Etot_1 += dE_1 + Vext(r_up, r_down, R, rho);

    dE_2 = 0.5 * (localEnergy1(Aup, Adown, Xold_up, Xold_down, r_up, r_down, eta[ee]) + localEnergy2(Aup, Adown, Xold_up, Xold_down, r_up, r_down, eta[ee]));
    Etot_2 += dE_2 + Vext(r_up, r_down, R, rho);
    // printf("mm: %d\t dE_1: %lf\t Etot_1: %lf\n",mm, dE_1, Etot_1 );
    }
    for (int jj = 0; jj < N; jj++){ // move all the 8 electrons at once
        for (int ii = 0; ii < DIM; ii++){
            Xnew[ii][jj] = Xold[ii][jj] + delta * (rand()/(double)RAND_MAX - 1.0/2.0);
        }
    }
    create_submat_and_r(Xnew, Xnew_up, Xnew_down, r_up, r_down);
    init_Slater_matrices(Aup, Adown, eta[ee], r_up, r_down, Xnew_up, Xnew_down);
    detAup = determinant_of_matrix(Aup); // Slater determinant
    detAdown = determinant_of_matrix(Adown); // Slater determinant
    init_Slater_matrices(Aup, Adown, eta[ee], r_up, r_down, Xnew_up, Xnew_down); // re-initialize
    psiNew = detAup * detAdown;
    
    w = pow(psiNew/psiOld, 2.0);
    random =  rand()/(double) RAND_MAX;
    // printf("w: %lf\t rand: %lf\n",w, random );
    
    if ( random <= w ){ // Metropolis
        /* Accept the move of the electron */
        for (int jj = 0; jj < N; jj++){
            for (int ii = 0; ii < DIM; ii++)
            {
               Xold[ii][jj] = Xnew[ii][jj];
            }
        }
    create_submat_and_r(Xold, Xold_up, Xold_down, r_up, r_down); 
    init_Slater_matrices(Aup, Adown, eta[ee], r_up, r_down, Xold_up, Xold_down); 
    // psiOld = psiNew;
    cont++;
    }
    
}// closes MC
printf("cont: %d\n", cont);
printf("Accept. prob: %lf\n", (double)cont/((M+thermal)));

dE_1 = localEnergy1(Aup, Adown, Xold_up, Xold_down, r_up, r_down, eta[ee]);
Etot_1 += dE_1 + Vext(r_up, r_down, R, rho);
dE_2 = 0.5 * (localEnergy1(Aup, Adown, Xold_up, Xold_down, r_up, r_down, eta[ee]) + localEnergy2(Aup, Adown, Xold_up, Xold_down, r_up, r_down, eta[ee]));
Etot_2 += dE_2 + Vext(r_up, r_down, R, rho);
Etot_1 /= M;
Etot_2 /= M;
Energies_1[ee] = Etot_1;
Energies_2[ee] = Etot_2;
printf(" Etot_1: %lf Etot_2: %lf\n", Etot_1, Etot_2 );
gsl_matrix_free(Aup);   gsl_matrix_free(Adown);
} // closes loop over eta

FILE *pf_t;
pf_t = fopen("dataEnEl.csv","w");
fprint_two_vec(pf_t, eta,  Energies_1, DIMe);
fclose(pf_t);



}