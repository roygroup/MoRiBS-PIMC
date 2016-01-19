#include "mc_utils.h"
#include "mc_confg.h"

#include <cmath>
#include <assert.h>
#include <stdlib.h>

//------------ memory allocation ----------------

double **doubleMatrix(long row, long col)
// based on nrutils.c dmatrix 
{
   double **mat_;

   mat_    = new double*[row];
   mat_[0] = new double [row*col]; 

   assert(mat_    != NULL);
   assert(mat_[0] != NULL);
	   
   for (long i=1; i<row; i++)
         mat_[i] = mat_[i-1]+col;

   return mat_;
}  

int **intMatrix(long row, long col)
// based on nrutils.c imatrix 
{
   int **mat_;

   mat_    = new int*[row];
   mat_[0] = new int [row*col]; 

   assert(mat_    != NULL);
   assert(mat_[0] != NULL);
   
	   
   for (long i=1; i<row; i++)
         mat_[i] = mat_[i-1]+col;

   return mat_;
}  
 
long int **long_intMatrix(long row, long col)
// based on nrutils.c imatrix 
{
   long int **mat_;

   mat_    = new long int*[row];
   mat_[0] = new long int [row*col]; 

   assert(mat_    != NULL);
   assert(mat_[0] != NULL);
   
	   
   for (long i=1; i<row; i++)
         mat_[i] = mat_[i-1]+col;

   return mat_;
}   

void free_doubleMatrix (double **mat) 
{
  delete [] mat[0];
  delete [] mat;
}

void free_intMatrix (int **mat) 
{
  delete [] mat[0];
  delete [] mat;
}

void free_long_intMatrix (long int **mat) 
{
  delete [] mat[0];
  delete [] mat;
}


//-----------------------------------------------

void nrerror(const char proc[],const char error_text[]) //  error handler 
// use cerr ?
{
    cout << endl;
    cout << "run-time error..." << endl;
    cout << proc<<":  "<<error_text << endl;
    cout << "...now exiting to system..."<< endl;
    cout << endl;
    exit(1);
}

//--------------cubic spline from Numerical recipes ----------

void spline(double *x, double *y, int n, double yp1, double ypn, double *y2)
// INPUT x[0..n-1]    - GRID (monotonically increasing)
//       y[0..n-1]    - tabulated function
//       
//       yp1 and ypn  - first derivatives of the interplotating function
//                      at 0 and (n-1)    
//
// OUTPUT       
//       y2 - second derivatives of the interpolating function
{
     double p,qn,sig,un, *u;

     u=new double [n];
     if (yp1>0.99e30)
	 y2[0]=u[0]=.0;
     else{
         y2[0]=-.5;
	 u[0]=(3./(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-yp1);
     
     }

     for (int i=1;i<(n-1);i++){
         sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
	 p=sig*y2[i-1]+2.;
	 y2[i]=(sig-1.)/p;
	 u[i]=(y[i+1]-y[i])/(x[i+1]-x[i])-(y[i]-y[i-1])/(x[i]-x[i-1]);
	 u[i]=(6.*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
     }

     if (ypn> 0.99e30)
	 qn=un=0.;
     else{
         qn=.5;
	 un=(3./(x[n-1]-x[n-2]))*(ypn-(y[n-1]-y[n-2])/(x[n-1]-x[n-2]));
     }

     y2[n-1]=(un-qn*u[n-2])/(qn*y2[n-2]+1.);

     for (int k=n-2;k>=0;k--)
	 y2[k]= y2[k]*y2[k+1]+u[k];

     delete [] u; 
}

void splint(double *xa, double *ya, double *y2a, int n, double x, double &y)
// INPUT xa[0..n-1]  - GRID (monotonically increasing)
//       ya[0..n-1]  - tabulated function
//       
//       y2a[0..n-1] - array of second derivatives from the spline function
//
//	 x - coordinate to interpolate function at
// OUTPUT       
//       y - cubic-spline interpolation 
// WARNING - no range checking       
{
    int klo,khi,k;
    double h,b,a;

    klo=0;                 // bisection method to find right place in table
    khi=n-1;
    while (khi-klo > 1) {
	k=(khi+klo)>>1;
	if (xa[k]>x) khi=k;
	else klo=k;
    }

    h=xa[khi]-xa[klo];
//  if (h== 0.) nrerror ("splint: Bad input ... ");
        
    if (khi== klo) nrerror ("splint:"," Bad input ... ");

    a=(xa[khi]-x)/h;
    b=(x-xa[klo])/h;
    y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.;
}

void init_spline(double *grid,double *data,double *sdata,int maxsize)
{	
   const char *_proc_=__func__;   // "init_spline";

   if (maxsize<2) nrerror(_proc_,ERR_INDEX_EXCEED);

   double drl = (grid[1]-grid[0]);
   double dpl = (data[1]-data[0])/drl;                  // first derivative (left)
   
   double drr = (grid[maxsize-1]-grid[maxsize-2]);
   double dpr = (data[maxsize-1]-data[maxsize-2])/drr;  // first derivative (right)

   spline(grid,data,maxsize,dpl,dpr,sdata);             
}


//-----------------------------------------------------------------------------

void  mmsort(double *dist, int *labels,int count) // sort dist and rearrange labels
// dist   [1...count]
// labels [1...count]
{
  int i,j;
  double dtmp;
  int    itmp;

  for (j=2;j<=count;j++) //j<=count - this is correct
  { 
     dtmp = dist[j];
     itmp = labels[j];

     i = j-1;

     while ((i>0) && (dist[i]>dtmp))
     {
         dist[i+1]   = dist[i];
	 labels[i+1] = labels[i];
	 i--;
     }

     dist[i+1]   = dtmp;
     labels[i+1] = itmp;
  }
}

double dmax(double x, double y)
{
    if (x >= y)
    {
        return x;
    }
    return y;
}

double dmin(double x, double y)
{
    if (x <= y)
    {
        return x;
    }
    return y;
}

