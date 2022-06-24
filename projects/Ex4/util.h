#ifndef __UTIL_H__
#define __UTIL_H__

#define DIMx 100
#define DIMe 5
#define N2 4 // N = 8 electrons (4 spin up and 4 spin down)
#define N 8
#define DIM 3
#define M 1e4// number of MC steps

static size_t size = 4;


/* ======================= FUNCTION HEADERS ======================= */
double Psi0(double eta, double r);
double Psi1(double eta, double r, double x);
void initialize_mat_contents(gsl_matrix *matrix, double eta, double r[4], double R[3][4]);
gsl_matrix *invert_a_matrix(gsl_matrix *matrix);
void print_mat_contents(gsl_matrix *matrix);
double determinant_of_matrix(gsl_matrix *matrix);
double Vext(double r, double R, double rho);
int randomGenerator(int low, int up);

void gradient_of_matrix(gsl_matrix *matrix, double eta, double r[4], double R[3][4], int component);
void laplacian_of_matrix(gsl_matrix *matrix, double eta, double r[4], double R[3][4] );
// double localEnergy(double psiT);
void GDtoDR_old(gsl_matrix *m1, gsl_matrix *m2, gsl_matrix *m3 , gsl_matrix *detInv, double GDratio[DIM][N2]);
void GDtoDR_new(gsl_matrix *A, gsl_matrix *m1, gsl_matrix *m2, gsl_matrix *m3 , gsl_matrix *detInv, double GDratio[DIM][N2]);
void LDtoDR(gsl_matrix *mm, gsl_matrix *detInv, double LDratio[N2]);
double localEnergy1( double LDtoDRup[N2], double LDtoDRdown[N2]);
double localEnergy2(double LDtoDRup[N2], double LDtoDRdown[N2], double GDtoDRup[DIM][N2], double GDtoDRdown[DIM][N2]);

void initialization_of_wf(gsl_matrix *m, gsl_matrix *gm_x, gsl_matrix *gm_y, gsl_matrix *gm_z, gsl_matrix *lm, double GDrat[DIM][N2] , double LDrat[4], double detm, double eta, double r[4], double R[3][4]);
void evolution_of_wf(gsl_matrix *m, gsl_matrix *inv_m, gsl_matrix *gm_x, gsl_matrix *gm_y, gsl_matrix *gm_z, gsl_matrix *lm, double GDrat[DIM][N2] , double LDrat[4], double detm, double eta, double r[4], double R[3][4]);


#endif