#ifndef _MC_ESTIM_H
#define _MC_ESTIM_H 1

#include <fstream>

void   InitMCEstims(void);
void   DoneMCEstims(void);
void  ResetMCEstims(void);

void SaveDensities1D(const char [], double);

void IOxyzAng(int,const char []);

void SaveDensities2D(const char [], double,int); //added by Hui Li
void SaveGraSum(const char fname [], double acount); //added by Toby Zeng
void SaveRho1D(const char fname [], double acount, int mode); //added by Toby Zeng
void SaveDensities3D(const char [], double,int); //added by Toby Zeng
void SaveRhoThetaChi(const char [], double,int); //added by Toby Zeng

double GetPotEnergy_Diff(void); //added by Hui Li

double GetPotEnergy_Densities(void);
double GetPotEnergy(void);
double GetKinEnergy(void);

double GetConfPoten_Densities(void); // HA test

// super densities

void GetExchangeLength(void);
void GetAreaEstimators(void);
void GetAreaEstim3D(int);

extern double _areas3DMFF[6];  // block accumulated area tensor estimator in dopant-fixed frame
extern double _inert3DMFF[9]; // block accumulated classical moment of inertia in dopant-fixed frame

extern double _areas3DSFF[6];  // block accumulated area tensor estimator in dopant-fixed frame
extern double _inert3DSFF[9]; // block accumulated classical moment of inertia in dopant-fixed frame

void SaveExchangeLength (const char [], double, long int);
void SaveAreaEstimators (const char [], double, long int);
void SaveAreaEstim3D (const char [], double, long int,int);

// permutation sampling
void GetPermutation(void);

// Rotations

double GetRotEnergy(void);
double GetRotE3D(void); // get rotational energy for nonlinear rotor, added by toby
double GetRotE3Dstep(int, int); // real step in loop of GetRotE3D, added by toby
double GetExactRotEnergy(int);

void GetRCF(void);
void SaveRCF(const char [],double,int);

void GetExactRCF(double *,int,int);

extern int  * _pflags;

extern int PrintXYZprl; // for printing instantaneous xyz coordinates and permutation table for each bloc

extern double ErotSQ; // asymmetric top rotational energy square estimator
extern double Erot_termSQ; // sum of square of each bead's rotational energy estimator

// MPI rotational variables
extern double srotchunk; // chunk summation of rotational energy
extern double srotsum; // total summation of rotational energy
#endif  // mc_estim.h
