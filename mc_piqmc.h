#ifndef _MC_PIQMC_H
#define _MC_PIQMC_H 1

void MCMolecularMove(int);
void MCBisectionMove(int,int);
void MCRotationsMove(int);
// Toby adds rotation move for nonlinear rotor
void MCRotations3D(int);
void MCRot3Dstep(int, int, int, int, double,double,double,double,double,int, int, double &, double &);
void MCRotLinStep(int,int,int,int,double,double,double,double,double &,double &);
void Reflect_MF_XZ(void);
void Reflect_MF_YZ(void);
void Reflect_MF_XY(void);
void RotSymConfig(void);

void MCMolecularMoveExchange(int);
void MCBisectionMoveExchange(int,int);

double PotEnergy(int,double **);
double PotEnergy(int,double **,int);

double PotRotEnergy(int,double **,int it);   
double PotRotE3D(int,double *,int it);   

extern double  **MCTotal;  // MC counters (total number of moves)
extern double  **MCAccep;  // MC counters (number of accepted moves)

// counters for parallel 3d rotation move.  They are supposed to be declared for all CPUs
extern double MCRotChunkAcp;    // total number of rotational moves for one chunk loop
extern double MCRotChunkTot;    // total accept number of rotational moves for one chunk loop
extern double MCRotTot; // the sum of all MCRotChunkTot from all CPUs
extern double MCRotAcp; // the sum of all MCRotChunkAcp from all CPUs


void ResetMCCounts(void);
void MemAllocMCCounts(void);
void MFreeMCCounts(void);

extern int PrintYrfl; // integer flag for printing reflected coordinates
extern int PrintXrfl; // integer flag for printing reflected coordinates
extern int PrintZrfl; // integer flag for printing reflected coordinates

#endif  //mc_pimc.h
