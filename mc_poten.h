#ifndef _MC_POTEN_H
#define _MC_POTEN_H 1

#include "mc_confg.h"

void     InitPotentials(void);
void     InitRotDensity(void);

void     DonePotentials(void);
void     DoneRotDensity(void);

double SRotDens(double gamma,int atype);       // rotational density matrix 
double SRotDensDeriv(double gamma ,int atype); // the derivs of the rotational density matrix 
double SRotDensEsqrt(double gamma ,int atype); // the 2nd derivs of the rotational density matrix 

double   SPot1D(double,int);        // spline interpolation for 1D potentials
double   LPot2D(double,double,int); // linear interpolation for 2D potentials
double   DLPot2D(double,double,int); // linear interpolation for 2D potential differences added by Hui Li
 
// ---------  MODELS OF INTERACTION ---------------

const int  PRIMITIVE = 0; 
const int  EFFECTIVE = 1; 

const char PMODEL[2][MAX_STRING_LENGTH] = {"PRIMITIVE",
                                           "EFFECTIVE"};

const char EXT_DIFF[] = ".diff";  // data file, potential differences added by Hui Li 
const char EXT_PRIM[] = ".pot";  // data file, primitive approximation
const char EXT_EFFC[] = ".eff";  // data file, supperpos approximation 
const char EXT_ROTD[] = ".rot";  // data file, rotational density matrix 
const char EXT_RHO[] = ".rho";  // data file, asymmetric top rotational density matrix 
const char EXT_ESQ[] = ".esq";  // data file, asymmetric top rotational energy square estimator

extern double * vtable;          // the table that contains 3D potential, added by Toby
extern int      Rgrd;  // # of radial grid points to store 3D potential
extern int      THgrd; // # of theta grid points to store 3D potential. It should be 181 at this moment.
extern int      CHgrd; // # of chi grid points to store 3D potential. It should be 91, 181, or 361 at this moment.
extern double   Rvmax; // maximum radius for 3D potential extrapolation
extern double   Rvmin; // minimum radius for 3D potential extrapolation
extern double   Rvstep; // radial increment for 3D potential extrapolation

//--------------------------------------------------
#endif  // mc_poten.h
