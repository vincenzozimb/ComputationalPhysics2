#ifndef __UTIL_H__
#define __UTIL_H__


#define DIMe 20
#define N2 4 // N = 8 electrons (4 spin up and 4 spin down)
#define N 8
#define DIM 3
#define M 80000 // number of MC steps
#define thermal 10000 // thermalisation of the system : make the system evolve without calculating the local energy

static size_t size = 4;


/* ======================= FUNCTION HEADERS ======================= */
double Vext(double rUp[N2], double rDown[N2], double R, double rho);

 void print_matN(double mat[DIM][N]);
 void print_matN2(double mat[DIM][N2]);
 void print_vecN2(double vec[N2]);

// int randomGenerator_int(int low, int up);
// double randomGenerator_double(double low, double up);

void init_positions(double mat[DIM][N], double delta);
void init_Slater_matrices(gsl_matrix *m1, gsl_matrix *m2, double eta, double r1[N2], double r2[N2], double R1[DIM][N2], double R2[DIM][N2]); 
void create_submat_and_r(double mat[DIM][N], double mat_up[DIM][N2], double mat_down[DIM][N2], double r_up[N2], double r_down[N2]);
double Psi0(double eta, double r);
double Psi1(double eta, double r, double x);
double orbital(double R[DIM][N2], double r[N2], double eta, int ii, int jj);
double grad_orbital(double R[DIM][N2], double r[N2], double eta, int ii, int jj, int kk);
void initialize_grad_mat(gsl_matrix* grad_mat[], double R[DIM][N2], double r[N2], double eta);
double laplacian_orbital(double R[DIM][N2], double r[N2] ,int ii, int jj, double eta);
void initialize_lap_mat(gsl_matrix* mat, double R[DIM][N2], double r[N2], double eta);
double localEnergy1(gsl_matrix* Aup, gsl_matrix* Adown, double Rup[DIM][N2], double Rdown[DIM][N2], double r1[N2], double r2[N2], double eta);
double localEnergy2(gsl_matrix* Aup, gsl_matrix* Adown, double Rup[DIM][N2], double Rdown[DIM][N2], double r1[N2], double r2[N2], double eta);

void initialize_mat_contents(gsl_matrix *matrix, double eta, double r[N2], double R[DIM][N2]);
gsl_matrix *invert_a_matrix(gsl_matrix *matrix);
void print_mat_contents(gsl_matrix *matrix);
double determinant_of_matrix(gsl_matrix *matrix);

void MC_integration();

// void gradient_of_matrix_x(gsl_matrix *matrix, double eta, double r[N2], double R[DIM][N2]);
// void gradient_of_matrix_y(gsl_matrix *matrix, double eta, double r[N2], double R[DIM][N2]);
// void gradient_of_matrix_z(gsl_matrix *matrix, double eta, double r[N2], double R[DIM][N2]);

// void laplacian_of_matrix(gsl_matrix *matrix, double eta, double r[N2], double R[DIM][N2] );

// void GDtoDR(gsl_matrix *m1, gsl_matrix *m2, gsl_matrix *m3 , gsl_matrix *detInv, double GDratio[DIM][N2]);
// void LDtoDR(gsl_matrix *mm, gsl_matrix *detInv, double LDratio[N2]);

// double localEnergy2(double LDtoDRup[N2], double LDtoDRdown[N2], double GDtoDRup[DIM][N2], double GDtoDRdown[DIM][N2]);

//void getRatio(gsl_matrix *matrix, gsl_matrix *detInv, double ratio[N2]);

#endif