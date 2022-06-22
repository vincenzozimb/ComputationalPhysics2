#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <string.h>

#include "print_routines.h"
#define DIMx 100
#define DIMe 5
#define N2 4 // N = 8 electrons (4 spin up and 4 spin down)
#define N 8
#define DIM 3
#define M 100000 // number of MC steps

static size_t size = 4;

/* ======================= FUNCTION HEADERS ======================= */
double Psi0(double eta, double r);
double Psi1(double eta, double r, double x);
void initialize_mat_contents(gsl_matrix *matrix, double eta, double r[4], double R[3][4]);
gsl_matrix * invert_a_matrix(gsl_matrix *matrix);
void print_mat_contents(gsl_matrix *matrix);
double determinant_of_matrix(gsl_matrix *matrix);
double Vext(double r, double R, double rho);

double localEnergy();

/* ======================= MAIN ======================= */
int main(){

/* Parameters of the system */
char atom[] = "Na"; // specify the atom type. Choose "Na" for sodium or "K" for potassium
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

for (int ii = 0; ii < DIMx; ii++)
{
    x[ii] = ii * h;
    Vex[ii] = Vext(x[ii],Ray, rho);
}

FILE *pf;
pf = fopen("dataVext.csv", "w");
fprint_two_vec(pf, x, Vex, DIMx);
fclose(pf);

double r2[4], r[4];
for (int jj = 0; jj < 4; jj++)
{
    r2[jj] = Xold[0][jj] * Xold[0][jj] + Xold[1][jj] * Xold[1][jj] + Xold[2][jj] * Xold[2][jj]; 
    r[jj] = sqrt(r2[jj]);
    //printf("r2[%d] = %lf\n", jj, r2[jj]);
}


/* Variational parameters */
double eta[DIMe];
double eta_init = 1.0;
double eta_end = 5.0;
double deta = (eta_end - eta_init)/(DIMe);

/* Energy and variance */
double Etot = 0.0;
double Etot2 = 0.0;

double gradAup[N2][N2];

/* Slater matrix */
gsl_matrix *Aup = gsl_matrix_alloc(size,size);
gsl_matrix *Adown = gsl_matrix_alloc(size,size);
gsl_matrix *Aup_inv = gsl_matrix_alloc(size,size);
gsl_matrix *Adown_inv = gsl_matrix_alloc(size,size);

double detAup, detAdown;

double psiOld, psiNew;

/* Evolution of the system : need 4 electrons to evolve individually */
/* Non interacting electrons - no Jastrow factor */
for (int ee = 0; ee < DIMe; ee++)
{
    eta[ee] = eta_init + ee * deta;

    Etot = 0.0;
    Etot2 = 0.0;

    /* Initialization of the matrix of positions : column -> electron, row -> x,y,z (determine th eposition in which
     I'm considering the electron )*/
    printf("Xold: \n");
    for (int ii = 0; ii < DIM; ii++)
    {
        for (int jj = 0; jj < N2; jj++)
        {
            Xold[ii][jj] = h * (rand()/(double)RAND_MAX - 1.0/2.0);
            printf("%lf\t", Xold[ii][jj]);
        }
        printf("\n");
    }
    printf("\n");
    /* Calculation of r */
    for (int ii = 0; ii < N2; ii++)
    {
        r2[ii] = Xold[0][ii] * Xold[0][ii] + Xold[1][ii] * Xold[1][ii] + Xold[2][ii] * Xold[2][ii]; 
        r[ii] = sqrt(r2[ii]);
    }
    
    /* Initialization of matrices A and A^(-1) */
    initialize_mat_contents(Aup, eta[ee], r, Xold);
    printf("Aup:\n");
    print_mat_contents(Aup);
    printf("\n");
    initialize_mat_contents(Adown, eta[ee], r, Xold);
    Aup_inv = invert_a_matrix(Aup);
    Adown_inv = invert_a_matrix(Adown);
     /* Re-initialization */
    initialize_mat_contents(Aup, eta[ee], r, Xold);
    initialize_mat_contents(Adown, eta[ee], r, Xold);
     /* Slater Determinant*/
    detAup = determinant_of_matrix(Aup);
    printf("det : %lf\n", detAup);
    detAdown = determinant_of_matrix(Adown);

    psiOld = detAup * detAdown;

    /* Monte Carlo loop */
    for (int mm = 0; mm < M; mm++)
    {
        // Sposto su un indice preso a random 
    
        //for (int ii = 0; ii < N2; ii++) // here am I making evolve the whole system?
        // {
            for (int jj = 0; jj < DIM; jj++)
            {
                Xnew[ii][jj] = Xold[ii][jj] + h * (rand()/(double)RAND_MAX - 1.0/2.0);
            }
            
        // }
        initialize_mat_contents(Aup, eta[ee], r, Xnew);
        initialize_mat_contents(Adown, eta[ee], r, Xnew);
        detAup = determinant_of_matrix(Aup);
        detAdown = determinant_of_matrix(Adown);
        // printf("Det: %lf\n", detAup);
        psiNew = detAup * detAdown;

        /* Metropolis test */
        if (rand()/(double)RAND_MAX <= psiNew * psiNew/(psiOld * psiOld))
        {
            for (int ii = 0; ii < N2; ii++)
            {
                for (int jj = 0; jj < DIM; jj++)
                {
                    Xold[ii][jj] = Xnew[ii][jj];
                }
            }
            psiOld = psiNew;
        }
        /* Compute local energy */
        // avoided thermalisation

    
    // Implementation of the local energy : laplacian and gradient
    }
    

}

// gsl_matrix_free(Aup);
// gsl_matrix_free(Adown);
// gsl_matrix_free(Aup_inv);
// gsl_matrix_free(Adown_inv);






}
/* ======================= FUNCTION BODIES ======================= */

double Psi0(double eta, double r){
    return exp(-r * r /(2.0 * eta * eta));
}

double Psi1(double eta, double r, double x){
    return x * exp(-r * r /(2.0 * eta * eta));
}

void initialize_mat_contents(gsl_matrix *matrix, double eta, double r[4], double R[3][4]){

    double val = 0.0; 
for ( size_t jj = 0; jj < size; jj++) // run over collumns the wf of the 4 electrons
    {
    for ( size_t ii = 0; ii < size; ii++) // run over row 4 positions of electrons
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

double detMat;
    
    gsl_permutation *p = gsl_permutation_alloc(size);
    int s;

    // Compute the LU decomposition of this matrix (s = 1)
     gsl_linalg_LU_decomp(matrix, p, &s);
     detMat = gsl_linalg_LU_det(matrix, s);

    gsl_permutation_free(p);
return detMat;
}

double Vext(double r, double R, double rho){
    double a = 2.0 * M_PI * rho;
    if (r <= R)
    {
        return a * (r * r/3.0 - R * R);
    } else {
        return a * -2.0/3.0 * R * R * R/r;
    
    }
    
}

double localEnergy(){
    
}