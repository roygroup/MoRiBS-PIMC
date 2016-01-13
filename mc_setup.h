#ifndef _MM_SETUP_H
#define _MM_SETUP_H 1

#include "mc_confg.h"
#include <string>
 
//------------- MC FLAGS ---------------------------------

extern bool    WORM;              // use the worm algorithm 

//------------- MC SYSTEM ---------------------------------

extern bool    IMPURITY;     // set true if there is a molecule in the system
extern bool    MINIMAGE;     // set true to apply minimum image convention 

extern bool    BOSONS;       // true if there're bosons in the system
extern int     BSTYPE;       // atom type for bosons

extern int     IMTYPE;       // atom type for dopant molecule

extern int     ISPHER;        // whether to treat the asymmetric top dopant as spherical particle. 0: no; 1: yes
extern int     IREFLY;        // whether to reflect all particles wrt the xz plane of the dopant, 0: no; 1: yes
extern int     IREFLX;        // whether to reflect all particles wrt the yz plane of the dopant, 0: no; 1: yes
extern int     IREFLZ;        // whether to reflect all particles wrt the xy plane of the dopant, 0: no; 1: yes
extern int     IROTSYM;        // whether to rotate the dopants by their body-fixed axis, 0: no; 1: yes
extern int     NFOLD_ROT;           // foldness of rotational symmetry of the dopant

extern bool    ROTATION;     // set to 1 to account for the rotational degrees of freedom

extern bool    FERMIONS;     // true if there're fermions in the system
extern int     FERMTYPE;     // atom type for fermions

extern int     NUMB_ATOMS;   // total number of atoms  
extern int     NUMB_MOLCS;   // total number of molecules  

extern int     NUMB_ATOMTYPES; // total number of atoms types 
extern int     NUMB_MOLCTYPES; // total number of molecules types

extern int     NDIM;
extern double  Temperature;
extern double  Density;
extern double  BoxSize;

extern double  MCBeta;     // imaginary time step 
extern double  MCTau;      // imaginary time step 

extern double  MCRotTau;   // imaginary time step for rotational degrees of freedom

extern int RotDenType; // Rotational Density Type, exact (default 0) or rattle and shake(1)
extern int RotOdEvn; // rotor symmetry info: -1 distinguishable, 0 pH2O type, 1 oH2O type
extern double RotEoff; // offset for the rotational energy estimator taken from Noya's formula at relative Euler angles = 0 0 0
extern double X_Rot; // rotational constant along principal x-axis in the unit of cm-1
extern double Y_Rot; // rotational constant along principal y-axis in the unit of cm-1
extern double Z_Rot; // rotational constant along principal z-axis in the unit of cm-1
extern int    RNratio; // ratio between RS and Noya steps in hybrid rotational energy estimation
extern int    NumbRotLim; // limit of number of one type of rotors

extern int NumbAtoms; // total number of atoms and molecules
extern int NumbTypes; // Number of particles' types

typedef struct TParticle
{
   int    numb;                    // number of atoms/molecule of this type  
   double mass;                    // mass of atom/molecule
   double brot;                    // rotational constant (for molecules only)  

   char   type[MAX_STRING_LENGTH]; // type of atoms/molecules
   char   fpot[MAX_STRING_LENGTH]; // the file with the tabulated potential
 
   int    stat;                    // statistics: Boltzmann,BE
   int    pmod;                    // model of interaction (1D/2D, primitive/effective)

   int    molecule; // 1 molecule, 0 - atoms
   double rtstep;   // step for the rotational degrees of freedom
// char   rdens [MAX_STRING_LENGTH];// the file name with the rotational density matrix

   int    offset;   // MCCoords[][offset+NumbTimes*atom+it]
   int    gatom;    // global atom counter: offset = NumbTimes*gatom

   double mcstep;   // MC step for molecular move
   int    levels;   // number of levels for multilevel Metropolis
   int    mlsegm;   // mlsegm = (int)pow(2.0,mclevels)

   double lambda;   // \hbar^2/2m
   double twave2;   // 4\lamda*tau
};

extern TParticle MCAtom[];

// -------------- MC TABLES -------------------------------   

extern int * MCType;   // convert atom number into atom type
extern int * PIndex;   // permutation index
extern int * RIndex;   // inverse permutation index

//------------- MC PARAMETERS -----------------------------

extern int NumbTimes;    // number of time slices (Max Number of Slices)
extern int MaxnTimes;    // Max number of time slices (NumbAtoms*NumbTimes)

extern int NumbRotTimes; // number of time slices for rotational degrees of freedom
extern int RotRatio;     // RotRatio = NumbTimes/NumbRotTimes;

extern long int NumberOfMCPasses;  // number of steps within a block
extern long int NumberOfMCBlocks;  // number of blocks
extern long int InitialBlock; // starting block number
extern long int NumberOfEQBlocks;  // number of equilibr blocks 

//------------- NONLINEAR MOLECULAR PARAMETERS--------------

const int SizeRotDen=181*361*361;
const int SizePotTab=501*181*181;

//------------- MPI PARAMETERS ----------------------------
extern int NProcs; // the number of processors as a global variable
extern int chunksize;  // the size of a chunk of rotational time slices treated by MPI
const int tagrho = 1;  // tag for passing rhoprp
const int tagero = 2;  // tag for passing erotpr
const int tagpot = 3;  // tag for passing potential
const int tagWORM = 4; // tag for passing WORM
const int tagMCType = 5; // tag for passing MCType[NumbAtoms]
const int tagMCAtom_molecule = 6; // tag for passing MCAtom[MAX_NUMBER_TYPES].molecule
const int tagRotRatio = 7; // tag for passing RotRatio
const int tagWormtype = 8; // tag for passing Worm.type
const int tagNumbRotTimes = 9; // tag for passing NumbRotTimes
const int tagMCAtom_offset = 10; // tag for passing MCAtom[MAX_NUMBER_TYPES].offset
const int tagIMTYPE = 11; // tag for passing IMTYPE
const int tagMCAtom_numb = 12; // tag for passing MCAtom[MAX_NUMBER_TYPES].numb
const int tagNumbTypes = 13; // tag for passing NumbTypes
const int tagNumbAtoms = 14; // tag for passing NumbAtoms
const int tagMCTau = 15; // tag for passing MCTau
const int tagTAGRUN = 16; // tag for passing tagrunning
const int tagROTATION = 17; // tag for passing the flag of ROTATION
const int tagMCStartBlock = 18; // tag for passing MCStartBlock
const int tagNumberOfMCBlocks = 19; // tag for passing NumberOfMCBlocks
const int tagNumberOfMCPasses = 20; // tag for passing NumberOfMCPasses
//const int tagTZMAT = 21; // tag for passing TZMAT
const int tagRotDenType = 21;
const int tagRotOdEvn = 22;
const int tagRotEoff = 23;
const int tagARot = 24;
const int tagBRot = 25;
const int tagCRot = 26;
const int tagMCRotTau = 27;

extern int tagrunning; // tag for the running MPI tag, which runs with the loop but alway within 0-32767
const int tagUpper = 32767; // MPI tag upper limit
extern long int SEED; // random seed that depends on CPU id

//------------- OpenMP Parameters ----------------
extern int NThreads; // the number of threads as a global variable

// number of MC steps to skip ...

extern int MCSKIP_RATIO;     //  to save information regarding the accept ratio
extern int MCSKIP_TOTAL;     //  to save accumulated average
extern int MCSKIP_AVERG;     //  to evaluate averages

// MC move types
const int MCMAXMOVES = 3;   // Max number of different types of MC moves
const int MCMOLEC    = 0;   // "molecule"  move
const int MCMULTI    = 1;   //  multilevel move
const int MCROTAT    = 2;   //  rotational degrees of freedom

//-------------MC STATUS -----------------------------------

extern long int MCStartBlock;

//------------ MC DATA STORAGE -----------------------------

extern double ** MCCoords;   // translational degrees of freedom
extern double ** MCCosine;   // orientational cosines
extern double ** MCAngles;   // cost and phi

//------------ Initial MCCoords and MCAngles;
extern double * MCCooInit;   // store the read in MCCoords
extern double * MCAngInit;   // store the read in MCAngles

//extern double ** TZMAT; // a temporary matrix for testing data structure

extern double ** newcoords;  // buffer for new coordinates
extern int     * atom_list;  // buffer for atom labels

extern double *  rhoprp;     // rotatinal propagator for non-linear rotor
extern double *  erotpr;     // rotational energy estimator for non-linear rotor
extern double *  erotsq;     // rotational energy square estimator for non-linear rotor

extern int       InitMCCoords; // integer flag for read in MCCoords;

void MCMemAlloc(void); // allocate memory
void MCMemFree(void);  // free memory

//------------ MC SYSTEM OF UNITS --------------------------

typedef struct TSystemOfUnits
{
   double mass;        
   double energy;       
   double length;
   double temperature;
// double momentum;
// double velocity;
// double time;

   string slength;
   string senergy;
};

extern TSystemOfUnits Units;

//-------------------------
extern int   MPIsize;    // MPI
extern int   MPIrank;    // MPI
//-----------------------------

void MCInitParams(void);
void MCSetUnits(void);
void MCInit(void);
void MCConfigInit(void);

void MCSetUnits_HO_TEST(void);

#endif  //MM_setup.h
