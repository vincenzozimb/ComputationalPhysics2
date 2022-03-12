#ifndef __UTIL_H__
#define __UTIL_H__

/* typedef */
typedef struct Param{
    double xi;
    double E;
}Param;


/* functions headers */
double potential(double x);
double F(double x, void *p);
void initialize_pot(double x[], double V[], int dim);
void save_data(double x[], double V[], complex double psi[], int dim, FILE *outfile);


#endif
