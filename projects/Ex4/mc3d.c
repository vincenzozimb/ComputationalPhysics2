
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "print_routines.h"

#define DIMx 100
#define DIMa 50

#define M 5e6 // # of MC steps

/* ======================= FUNCTION HEADERS ======================= */

double PsiT(double alpha, double x);
double localEnergy(double alpha, double x);
double localEnergy2(double alpha, double x0, double h);
void PsiT_ex(double alpha, double x[DIMx], double psi[DIMx]);
double normCalculator(double psi[], double dx, int dim);


/* ======================= MAIN ======================= */

int main(){
    
/* Space parameters */
    double L = 20.0; // cube
    double h = L / (DIMx - 1.0);
    double posOld = -L/2.0;
    double posNew = -L/2.0;
    double x[DIMx];
    for (int ii = 0; ii < DIMx; ii++)
    {
           x[ii] = - L/2.0 + h * ii;
    }

    double X[DIMx][3];
    double r2[DIMx];
    for (int ii = 0; ii < DIMx; ii++)
    {
        for (int jj = 0; jj < 3; jj++)
        {
            X[ii][jj] = - L/2.0 + h * ii;
        }
        r2[ii] = X[ii][0] * X[ii][0] + X[ii][1] * X[ii][1] +X[ii][2] * X[ii][2];
    }
    
/* Variational parameters */
    double alpha_init = 0.3;
    double alpha_end = 3.0;

    double dalpha = (alpha_end - alpha_init)/(DIMa - 1.0);
    printf("dalpha = %lf\n", dalpha);   
    double alpha[DIMa];
    
/* Wave functions */
    double psiOld, psiNew;

/* Energy and variance*/
    double Etot, Etot2, dE;
    double Var;
    double Energies[DIMa];
    double Variances[DIMa];

    double ENERGY[DIMa];
    double VARIANCE[DIMa];

/* Implementation of MonteCarlo algorithm */
// Because of isotropy, work with 1D harmonic oscillator and at the end I multiply for 3
        for (int ii = 0; ii < DIMa; ii++){
           
        alpha[ii] = alpha_init + ii * dalpha;
        
        Etot = 0.0; 
        Etot2 = 0.0;

        posOld = posOld + h * (rand()/(double)RAND_MAX - 1.0/2.0); // add PosOld
        psiOld = PsiT(alpha[ii], posOld);

        /* Monte Calro loop */
        for (int mm = 0; mm < M; mm++)
        {
            posNew = posOld + h * (rand()/(double)RAND_MAX  - 1.0/2.0);
            psiNew = PsiT(alpha[ii], posNew);

            /* Metropolis */
            if (rand()/(double)RAND_MAX <= psiNew * psiNew/(psiOld * psiOld))
            {
                posOld = posNew;
                psiOld = psiNew;
            }   
             dE = localEnergy(alpha[ii], posOld); 
            // dE = localEnergy2(alpha[ii], posOld, h);
            Etot += dE;
            Etot2 += dE * dE;
        }
        /* Compute local energy */
        Etot /= M;
        Etot2 /= M;
        Var = Etot2 - Etot *Etot;
        Energies[ii] = Etot;
        Variances[ii] = Var;
    }


    for (int ii = 0; ii < DIMa; ii++)
    {  
        ENERGY[ii] = 3.0 * Energies[ii];
        VARIANCE[ii] = 3.0 * Variances[ii];
    }

/* Exact results from analytical calculations */
/* Exploit isotropy of ho */

double psiT[DIMx];
double norm[DIMa];
double El[DIMx][DIMa];
double p[DIMx][DIMa];
double Etot_ex[DIMa];
double Var_ex[DIMa];

    // for (int jj = 0; jj < DIMa; jj++)
    // {
    //     PsiT_ex(alpha[jj], x, psiT);
    //     norm[jj] = normCalculator(psiT, h, DIMx); 

    //     for (int ii = 0; ii < DIMx; ii++)
    //         {
    //             p[ii][jj] = pow(alpha[jj], 3.0) * exp(- alpha[jj] * alpha[jj] * r2[ii]) /(pow(M_PI, 3.0/2.0));
    //             El[ii][jj] = -1.0/2.0 * (-3.0 * alpha[jj] * alpha[jj] + r2[ii] * (pow(alpha[jj], 4.0) - 1.0));
    //         }
    // }

    // for (int jj = 0; jj < DIMa; jj++)
    // {   Etot_ex[jj] = 0.0;
    //     for (int ii = 0; ii < DIMx; ii++)
    //     {
    //         Etot_ex[jj] += p[ii][jj] * El[ii][jj]; // total energy
    //     }
    //     Etot_ex[jj] *= h;   
    // }

    for (int ii = 0; ii < DIMa; ii++)
    {
        Etot_ex[ii] = 3.0/4.0 * (alpha[ii]* alpha[ii] + 1.0 / (alpha[ii] * alpha[ii]));
    }

/* Print data on file */
FILE *pf, *pf1;
pf = fopen("dataMC3d.csv", "w");
pf1 = fopen("dataEx3d.csv", "w");

fprint_three_vec( pf,alpha, ENERGY, VARIANCE, DIMa );
fprint_three_vec(pf1, alpha, Etot_ex, Var_ex, DIMa);

fclose(pf);
fclose(pf1);
        
}

/* ======================= FUNCTION BODIES ======================= */
double PsiT(double alpha, double x){

    return exp(-1.0/2.0 * alpha * alpha * x * x);
    
}

double localEnergy(double alpha, double x){
    
     return  1.0/2.0 * (alpha * alpha + x * x *(1.0 - pow(alpha, 4.0)));

}

double localEnergy2(double alpha, double x0, double h){
    double df2;
    /* Second derivative */
    df2 = (PsiT(alpha, x0+h) + PsiT(alpha, x0-h) - 2.0 * PsiT(alpha, x0))/(h * h);
    return ( - 1.0/2.0 * df2  + 1.0/2.0 * x0 * x0 * PsiT(alpha, x0)) / PsiT(alpha, x0);

}

void PsiT_ex(double alpha, double x[DIMx], double psi[DIMx]){
    for (int ii = 0; ii < DIMx; ii++)
    {
         psi[ii] = exp(-1.0/2.0 * alpha * alpha * x[ii] * x[ii]);
    }
}

double normCalculator(double psi[], double dx, int dim){
    /* calculate norm with middle point integration*/
    double norm = 0.0;
    
    for(int i=0; i<dim; i++){
        norm += psi[i] * psi[i];
    }
    norm *= dx;
    norm = sqrt(norm);

    return norm * norm;
}

