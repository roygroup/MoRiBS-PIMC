#ifndef _MC_UTILS_H
#define _MC_UTILS_H 1

//------------ memory allocation ----------------

double      **doubleMatrix (long row, long col);
int         **intMatrix (long row, long col);
long int    **long_intMatrix (long row, long col);

void   free_doubleMatrix (double **mat);
void   free_intMatrix (int **mat); 
void   free_long_intMatrix (long int **mat); 

//-----------------------------------------------

const char ERR_INDEX_EXCEED[] =  "Index exceeds array size";

void nrerror(const char [],const char []); //  error handler 

//---- cubic spline from Numerical recipes ------

void spline(double *x, double *y, int n, double yp1, double ypn, double *y2);
void splint(double *xa, double *ya, double *y2a, int n, double x, double &y);

//------------------------------------------------

void init_spline(double *grid,double *data,double *sdata,int maxsize);

void  mmsort(double *, int *,int);  // sort dist and rearrange labels

double dmax(double,double);
double dmin(double,double);

#endif  // mc_utils.h
