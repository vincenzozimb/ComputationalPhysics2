#include <print_routines.h>

#include <stdio.h>
#include <assert.h>

void fprint_double(FILE *file, double d){
	fprintf(file, "%lf\t", d);
}

void fprint_double_newline (FILE *file, double d) {
	fprintf(file, "%lf\n", d);
}

void fprint_vec(FILE *file, double vec[], int dim){
	for(int i=0;i<dim;i++){
		fprintf(file,"%lf\n",vec[i]);
	}
	fprintf(file,"\n");
}

void fprint_mat(FILE *file, const double *A, int nrows, int ncols) {
	assert(nrows >= 0 && ncols >= 0);
	for (int i = 0; i < nrows; i++) {
		for (int j = 0; j < ncols; j++) {
			fprintf(file, "%lf\t", *(A + i * ncols + j));
		}
		fprintf(file, "\n");
	}
	fprintf(file, "\n");
}

void fprint_mat_complex(FILE *file, _Complex double (*A), int nrows, int ncols) {
	assert(nrows >= 0 && ncols >= 0);
	for (int i = 0; i < nrows; i++) {
		for (int j = 0; j < ncols; j++) {
			fprintf(file, "%.3lf+%.3lfi\t", creal(*(A + i * ncols + j)), cimag(*(A + i * ncols + j)));
		}
		fprintf(file, "\n");
	}
	fprintf(file, "\n");
}

void fprint_two_vec(FILE *file, double v1[], double v2[], int dim){
	assert(dim>0);
	for(int i=0;i<dim;i++){
		fprintf(file,"%lf\t",v1[i]);
		fprintf(file,"%lf\n",v2[i]);
	}
}