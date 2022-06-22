/* ----- Implementation of MC simulation for the calculation of the ground state energy of the 1D harmonic oscillator ----  */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "print_routines.h"

#define DIMx 100
#define DIMa 200

#define M 100000 // # of MC steps

/* ======================= FUNCTION HEADERS ======================= */

double PsiT(double alpha, double x);
double localEnergy(double alpha, double x);
void PsiT_ex(double alpha, double x[DIMx], double psi[DIMx]);
double normCalculator(double psi[], double dx, int dim);
int minimum();

/* ======================= MAIN ======================= */

int main(){

/*-------------- main parameters ---------------- */
    double L = 20.0;
    double posOld = -L/2.0;
    double posNew = -L/2.0;

    double h;
    h = L/(DIMx - 1.0);
    printf("h = %lf\n", h);

    double x[DIMx];
    for (int ii = 0; ii < DIMx; ii++)
    {
        x[ii] = -L/2.0 + h * ii;
    }

    double alpha_init = 0.5;
    double alpha_end = 1.5 ; 
    double alpha[DIMa];
    double dalpha;
    dalpha = (alpha_end - alpha_init)/(DIMa - 1.0);
    printf("dalpha = %lf\n", dalpha);

    double psiOld, psiNew;
    double psiT[DIMx];
    double norm[DIMa];

    double Etot = 0.0, Etot2 = 0.0, dE;
    double var, err;
    double Energies[DIMa], Variances[DIMa];

    double p[DIMx][DIMa];
    double El[DIMx][DIMa];
    double Etot_ex[DIMa];
    double Var_ex[DIMa];


/* Implementation of MonteCarlo algorithm */
    for (int ii = 0; ii < DIMa; ii++)
    {
        alpha[ii] = alpha_init + ii * dalpha;
        Etot = 0.0;
        Etot2 = 0.0;

        posOld = posOld + h * (rand()/(double)RAND_MAX - 1.0/2.0); // add posOld
        psiOld = PsiT(alpha[ii], posOld);

        /* MonteCarlo loop */
        for (int mm = 0; mm < M; mm++)
        {
            posNew = posOld + h * (rand()/(double)RAND_MAX  - 1.0/2.0);
            psiNew = PsiT(alpha[ii], posNew);

            /* Metrolopis to see if I accept the move or not */
            if (rand()/(double)RAND_MAX <= psiNew * psiNew/(psiOld * psiOld))
            {
                posOld = posNew;
                psiOld = psiNew;
            }   
        dE = localEnergy(alpha[ii], posOld);
        Etot += dE;
        Etot2 += dE * dE;
        }
        
        Etot /= M;
        Etot2 /= M;
        var = Etot2 - Etot *Etot;
        err = sqrt(var / M);
        Energies[ii] = Etot;
        Variances[ii] = var;
    }

/* Exact results from analytical caluclations */
    for (int jj = 0; jj < DIMa; jj++)
    { 
        PsiT_ex(alpha[jj], x, psiT);
        norm[jj] = normCalculator(psiT, h, DIMx); 

        for (int ii = 0; ii < DIMx; ii++)
        {
            p[ii][jj] = psiT[ii]*psiT[ii]/norm[jj]; // PDF
            El[ii][jj] = 1.0/2.0 * (alpha[jj] * alpha[jj] + x[ii] * x[ii] *(1.0 - pow(alpha[jj], 4.0))); // Local density
        }
        
    }
    
    // for (int jj = 0; jj < DIMa; jj++)
    // {   Etot_ex[jj] = 0.0;
    //     for (int ii = 0; ii < DIMx; ii++)
    //     {
    //         Etot_ex[jj] += p[ii][jj] * El[ii][jj]; // total energy
    //     }
    //     Etot_ex[jj] *= h;   

    //     Var_ex[jj] = 1.0/4.0 * (1.0 + pow((1.0 - pow(alpha[jj], 4.0)), 2.0) * 3.0/(4.0 * pow(alpha[jj], 4.0))) - Etot_ex[jj]* Etot_ex[jj];
    // }

    for (int ii = 0; ii < DIMa; ii++)
    {
        Etot_ex[ii] = 1.0/4.0 * (alpha[ii] * alpha[ii] + 1.0/(alpha[ii] * alpha[ii]));
        Var_ex[ii] = 1.0/4.0 * (1.0 + pow((1.0 - pow(alpha[ii], 4.0)), 2.0) * 3.0/(4.0 * pow(alpha[ii], 4.0))) - Etot_ex[ii]* Etot_ex[ii];
    }
    
    

/* Print data on files */
    FILE *pf, *pf1;
    pf = fopen("dataMC.csv", "w");
    pf1 = fopen("dataEx.csv", "w");
    
    fprint_three_vec(pf, alpha, Energies, Variances, DIMa);
    fprint_three_vec(pf1, alpha, Etot_ex, Var_ex, DIMa);
    
    fclose(pf);
    fclose(pf1);

/* Minimization algorithm to find the best alpha */
// Doesn't really work

    int pos_min;
    pos_min = minimum(Etot_ex);
    printf("alpha min : %lf\n", alpha[pos_min]);

}

/* ======================= FUNCTION BODIES ======================= */
double PsiT(double alpha, double x){

    return exp(-1.0/2.0 * alpha * alpha * x * x);
    
}

double localEnergy(double alpha, double x){
    
    return  1.0/2.0 * (alpha * alpha + x * x *(1.0 - pow(alpha, 4.0)));
    
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


int minimum(double E[DIMa]){

    int pos = 0;
    int ris = E[0];
    for (int ii = 1; ii < DIMa; ii++)
    {
       if (E[ii] != 0.0 && E[ii]< ris)
       {
            ris = E[ii];
            pos = ii; 
       }
       
    }
    return pos;
}
