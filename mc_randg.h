#ifndef _MC_RANDG_H
#define _MC_RANDG_H 1

double rnd (void);  // double random number generators in the range (0,1)
double rnd1(void); 
double rnd2(void); 
double rnd3(void); 
double rnd4(void); 
double rnd5(void); 
double rnd6(void); 
double rnd7(void); 
 
int nrnd1(int n);   // int random number generators in the range (0,n-1)
int nrnd2(int n);   
int nrnd3(int n);   
int nrnd4(int n);   

double gauss(double);

void RandomInit(int,int);          //  sprng init 
void RandomFree(void);             //  free memory used by sprng
void RandomIO(int, const char[]);  // check point for sprng streams

#endif  // mc_randg.h
