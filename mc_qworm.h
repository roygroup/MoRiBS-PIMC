#ifndef _MC_QWORM_H
#define _MC_QWORM_H 1

#include "mc_confg.h"

extern int oc_flag;    // 0 - open, 1 - close steps  DEBUG

extern double countCLOSE;
extern double countCLOSEA;
extern double countOPEN;
extern double WPOT_O;
extern double WKIN_O;

extern double WPOT_C;
extern double WKIN_C;

extern double WPOT_C_ACCEPT;
extern double WKIN_C_ACCEPT;

extern int rec;
extern int rec1;

extern int adv;
extern int adv1;

void MCWormInit(void);
void MCWormDone(void);
void MCWormMove(void);

bool WorldLine(int, int);            

typedef struct TPathWorm
{
   char    stype[MAX_STRING_LENGTH]; // type of atoms/molecules
   int     type;                     // atoms type to apply the worm algorithm (He)

   int     exists;  // set to 1 for G-space, to 0 for Z-space

   int     ira;
   int     masha;

   int     atom_i;  // ira's   world line
   int     atom_m;  // masha's world line    
 
   double  c;       // control relative stat of Z and G sectors    
   int     m;       // m tilde
};

extern TPathWorm Worm;

extern double **QWTotal;    // total number of worm moves
extern double **QWAccep;    // total number of accepted worm moves)

extern double countQW;      // total count of worm moves

// worm move types
const int QWMAXMOVES = 7;   // Max number of different types of worm moves

const int QW_OPEN    = 0;   
const int QW_CLOSE   = 1;   
const int QW_INSERT  = 2;   
const int QW_REMOVE  = 3;   
const int QW_ADVANCE = 4;   
const int QW_RECEDE  = 5;   
const int QW_SWAP    = 6; 
 
void ResetQWCounts(void);
void MemAllocQWCounts(void);
void MFreeQWCounts(void);
 
#endif  // mc_qworm.h
