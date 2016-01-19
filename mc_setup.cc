#include "mc_setup.h"
#include "mc_utils.h"
#include "mc_confg.h"
#include "mc_const.h"

#include <cmath>

//------------- MC FLAGS ---------------------------------

bool    WORM     = false;      // use the worm algorithm 

//------------- MC SYSTEM ---------------------------------

bool    IMPURITY = false;    // set true if there is a molecule in the system
bool    MINIMAGE = false;    // set true to apply minimum image convention 

int     IMTYPE   = -1;       // atom type for dopant molecule (rotational degrees of freedom)

int     ISPHER   = 0;        // whether to treat the asymmetric top dopant as spherical particle. 0: no; 1; yes
int     IREFLY   = 0;        // whether to reflect all particles wrt the xz plane of the dopant, 0: no; 1: yes
int     IREFLX   = 0;        // whether to reflect all particles wrt the yz plane of the dopant, 0: no; 1: yes
int     IREFLZ   = 0;        // whether to reflect all particles wrt the xy plane of the dopant, 0: no; 1: yes
int     IROTSYM  = 0;        // whether to rotate the dopants by their body-fixed axis, 0: no; 1: yes
int     NFOLD_ROT;           // foldness of rotational symmetry of the dopant

bool    ROTATION = false;    // set to 1 to account for the rotational degrees of freedom

bool    BOSONS   = false;    // true if there're bosons in the system
int     BSTYPE   = -1;       // atom type for bosons
 
bool    FERMIONS = false;    // true if there're fermions in the system
int     FERMTYPE = -1;       // atom type for fermions

int     NUMB_ATOMS = 0;      // total number of atoms  
int     NUMB_MOLCS = 0;      // total number of molecules  

int     NUMB_ATOMTYPES = 0;  // total number of atoms types 
int     NUMB_MOLCTYPES = 0;  // total number of molecules types
 
int     NDIM;
double  Temperature;
double  Density;
double  BoxSize;

double  MCBeta;
double  MCTau;      // imaginary time step 

double  MCRotTau;   // imaginary time step for rotational degrees of freedom

int RotDenType = 0; // Rotational Density Type, exact (default 0) or rattle and shake(1)
int RotOdEvn = 0; // rotor symmetry info: -1 distinguishable, 0 pH2O type, 1 oH2O type
double RotEoff = 0.0; // offset for the rotational energy estimator taken from Noya's formula at relative Euler angles = 0 0 0
double X_Rot = 0.0; // rotational constant A in the unit of cm-1
double Y_Rot = 0.0; // rotational constant B in the unit of cm-1
double Z_Rot = 0.0; // rotational constant C in the unit of cm-1
int    RNratio = 1; // ratio between RS and Noya steps in hybrid rotational energy estimation

int    NumbRotLim = 10; // limit of number of one type of rotors

int     NumbAtoms;  // total number of atoms and molecules
int     NumbTypes;  // Number of particles' types

TParticle MCAtom[MAX_NUMBER_TYPES];  // size should be NumbTypes+1

// -------------- MC TABLES -------------------------------   

int  * MCType;     // convert atom number into atom type
int  * PIndex;     // permutation index
int  * RIndex;     // inverse permutation index

//------------- MC PARAMETERS -----------------------------

int NumbTimes;     // number of time slices (number of slices per atom)
int MaxnTimes;     // NumbTimes*NumbAtoms - total number of slices 

int NumbRotTimes;  // number of time slices for rotational degrees of freedom
int RotRatio;      // RotRatio = NumbTimes/NumbRotTimes;

long int NumberOfMCPasses; // number of steps within a block
long int NumberOfMCBlocks; // number of blocks
long int InitialBlock; // starting block number
long int NumberOfEQBlocks; // number of equilibration blocks 

// number of MC steps to skip ...

int MCSKIP_RATIO = 100000;    //  to save information regarding the accept ratio
int MCSKIP_TOTAL = 10000;     //  to save accumulated average
int MCSKIP_AVERG = 1;         //  to evaluate averages

//------------- MC STATUS ----------------------------------

long int MCStartBlock;

//---------------- MPI PARALL ---------------
int NProcs; // the number of processors as a global variable
int chunksize;  // the size of a chunk of rotational time slices treated by MPI
int tagrunning; // tag for the running MPI tag, which runs with the loop but alway within 0-32767
long int SEED; // random seed that depends on CPU id

//---------------- OpenMP PARLL --------------
int NThreads; // the number of threads as a global variable

//------------ MC DATA STORAGE -----------------------------

double ** MCCoords;   // translational degrees of freedom
double ** MCCosine;   // orientational cosines
double ** MCAngles;

//------------ Initial MCCoords and MCAngles;
double * MCCooInit;   // store the read in MCCoords
double * MCAngInit;   // store the read in MCAngles

//------------ test the data structure of MCCoords ---------
//double ** TZMAT; // just a matrix

double ** newcoords;  // buffer for new coordinates
int    *  atom_list;  // buffer for atom labels

double *  rhoprp;     // rotatinal propagator for non-linear rotor
double *  erotpr;     // rotational energy estimator for non-linear rotor
double *  erotsq;     // rotational energy square estimator for non-linear rotor

int       InitMCCoords; // integer flag for read in MCCoords;

//----------------------------------------------------------

void initLattice_config(double **);
void replInitial_config(double **);

void MCMemAlloc(void)  // allocate  memmory 
{
   MCCoords  = doubleMatrix (NDIM,NumbAtoms*NumbTimes);  
   newcoords = doubleMatrix (NDIM,NumbAtoms*NumbTimes); 

// TZMAT = doubleMatrix (NDIM,NDIM);
 
   atom_list = new int [NumbAtoms];

// MCCosine  = doubleMatrix (NDIM,NumbAtoms*NumbTimes);  
   MCCosine  = doubleMatrix (3,NumbAtoms*NumbTimes); 
   MCAngles  = doubleMatrix (3,NumbAtoms*NumbTimes); 

   MCCooInit = new double [NDIM*NumbAtoms*NumbTimes];
   MCAngInit = new double [NDIM*NumbAtoms*NumbTimes];

// non-linear rotor
   rhoprp    = new double [SizeRotDen];
   erotpr    = new double [SizeRotDen];
   erotsq    = new double [SizeRotDen];

// TABLES
  
   MCType    = new int [NumbAtoms];
   PIndex    = new int [NumbAtoms];
   RIndex    = new int [NumbAtoms];
}

void MCMemFree(void)  //  free memory
{
   free_doubleMatrix(MCCoords);  
   free_doubleMatrix(newcoords); 

   delete [] atom_list;

   delete [] rhoprp;
   delete [] erotpr;
 
   free_doubleMatrix(MCCosine); 
   free_doubleMatrix(MCAngles); 

   delete MCType;
   delete PIndex;
   delete RIndex;
}

//------------ MC SYSTEM OF UNITS --------------------------

TSystemOfUnits Units;

//-------------------------
int   MPIsize;    // MPI
int   MPIrank;    // MPI
//-----------------------------

void MCSetUnits(void)
{
   Units.temperature = 1.0;   // Kelvin
   Units.energy      = 1.0;   // Kelvin
   Units.length      = 1.0;   // Angstrom
   Units.mass        = 1.0;   // amu

   Units.senergy     = "Kelvin";
   Units.slength     = "Angstrom";

   Temperature      /= Units.temperature;
   Density          *= (Units.length*Units.length*Units.length);

   double lambda     = 100.0*(HBAR*HBAR)/(AMU*K_B);   // \AA^2 K

   for (int atype=0;atype<NumbTypes;atype++)
   {
      MCAtom[atype].mcstep /= Units.length;
      MCAtom[atype].mass   /= Units.mass;
      MCAtom[atype].brot   /= Units.energy;

      MCAtom[atype].lambda  = 0.5*lambda/MCAtom[atype].mass;  // (hbar^2/2m)
   } 
}

void MCSetUnits_HO_TEST(void)
{
   Units.temperature = 1.0;   // Kelvin
   Units.energy      = 1.0;   // Kelvin
   Units.length      = 1.0;   // Angstrom
   Units.mass        = 1.0;   // amu

   Units.senergy     = "arb units";
   Units.slength     = "arb units";

   Temperature      /= Units.temperature;
   Density          *= (Units.length*Units.length*Units.length);

   double lambda     = 1.0;   // \AA^2 K

   for (int atype=0;atype<NumbTypes;atype++)
   {
      MCAtom[atype].mcstep /= Units.length;
      MCAtom[atype].mass   /= Units.mass;
      MCAtom[atype].brot   /= Units.energy;

      MCAtom[atype].lambda  = 0.5*lambda/MCAtom[atype].mass;  // (hbar^2/2m)
   } 
}

void MCInitParams(void)
{
   const char *_proc_=__func__;    // "MCInitParams()";

   double mass;
   double brot;

   for (int atype=0;atype<NumbTypes;atype++)
   {
       string stype=MCAtom[atype].type;

       if  (stype == HE4)
       {
          mass = MASS_HE4; 
          brot = 0.0; 
       }
       else
       if (stype ==H2)
       {
          mass=MASS_H2;
          brot = 0.0; 
       }
       else
       if  (stype == OCS)
       {
          mass =(MASS_O16 + MASS_C12 + MASS_S32);
          brot = B_OCS;
       } 
       else
       if  (stype == N2O)
       {
          mass = (2.0*MASS_N14 + MASS_O16);
          brot = B_N2O;
       } 
       else 
       if  (stype == CO2)
       {
          mass =(MASS_C12 + 2.0*MASS_O16);
          brot = B_CO2;
       } 
       else
       if  (stype == CO)
       {
          mass =(MASS_C12 + MASS_O16);
          brot = B_CO;
       }
       else
       if  (stype == HCN)
       {
          mass =(MASS_H1 + MASS_C12 + MASS_N14);
          brot = B_HCN;
       } 
       else 
       if  (stype == HCCCN)
       {
          mass =(MASS_H1 + 3.0*MASS_C12 + MASS_N14);
          brot = B_HCCCN;
       } 
       else
       if  (stype == H2O)
       {
          mass =(2.0*MASS_H1 + MASS_O16);
       }
       else
       if  (stype == SO2)
          mass =(2.0*MASS_O16 + MASS_S32);
       else
       if  (stype == HCOOCH3)
          mass =(4.0*MASS_H1 + 2.0*MASS_O16 + 2.0*MASS_C12);
       else 
       nrerror(_proc_,"Unknown atom/molecule type");

       MCAtom[atype].mass = mass;
       MCAtom[atype].brot = brot;
   }
}

void MCInit(void)  // only undimensional parameters in this function 
{
   const char *_proc_=__func__;    // "MCInit()";

//  INITIALIZE MC TABLES ---------------------------------------

   int natom = 0;    // map atom number into atom type
   for (int type=0;type<NumbTypes;type++)
   for (int atom=0;atom<MCAtom[type].numb;atom++)
   {
      MCType[natom] = type;
      natom ++;
   }

   for (int atom=0;atom<NumbAtoms;atom++)
   { 
      PIndex[atom] = atom;
      RIndex[atom] = atom;
   }
// ------------------------------------------------------------

// BoxSize  =  pow((double)NumbAtoms/Density,1.0/(double)NDIM); 
// define a box size based on number of atoms only (molecules excluded)
  
   BoxSize  =  pow((double)(NUMB_ATOMS+NUMB_MOLCS)/Density,1.0/(double)NDIM);

   MCBeta   =  1.0/Temperature;
   MCTau    =  MCBeta/(double)NumbTimes;

   if (ROTATION)
   MCRotTau =  MCBeta/(double)NumbRotTimes;

   RotRatio  = 1;  // div_t quot - it's important for the area estimator
                   // even without rotations

   if (ROTATION)
   {
      RotRatio = NumbTimes / NumbRotTimes;  // div_t quot
      int rt   = NumbTimes % NumbRotTimes;  // div_t rem
#ifndef ROTS_TEST
      if (rt)
      nrerror (_proc_,"NumbTimes is not proportional to NumbRotTimes");
#endif
   }

   for (int type=0;type<NumbTypes;type++)  
   {
      MCAtom[type].twave2 = 4.0*MCAtom[type].lambda * MCTau;   // thermal wavelength squared
      MCAtom[type].mlsegm = (int)pow(2.0,MCAtom[type].levels); // segmen size for multilevel
      
      if (MCAtom[type].mlsegm >= NumbTimes)
      nrerror (_proc_,"Segment size is larger then a number of time slices");
   } 

   int bcount = 0;  // number of bosons' types
   int fcount = 0;  // number of fermions' types
   int icount = 0;  // number of molecules (impurities)

   BOSONS   = false;
   FERMIONS = false;

   for (int type=0;type<NumbTypes;type++)
   {
      if (MCAtom[type].stat == BOSE)
      {
         BOSONS = true;
         BSTYPE = type;
         bcount ++;  
      } 

      if (MCAtom[type].stat == FERMI)
      {
         FERMIONS = true;
         FERMTYPE = type;
         fcount ++;  
      } 

      if ((MCAtom[type].molecule == 1)||(MCAtom[type].molecule == 2))
      {
         IMTYPE = type;
         icount ++;  
      } 
   }

   if (bcount>1)
   nrerror (_proc_,"Too many boson atoms' types");

   if (fcount>1)
   nrerror (_proc_,"Too many fermion atoms' types");

   if (icount>1)
   nrerror (_proc_,"Too many dopant molecule' types");

   if ((icount == 0) && (IMPURITY))
   nrerror (_proc_,"Impurity type not defined");

   if ((BOSONS && FERMIONS) && (BSTYPE == FERMTYPE))
   nrerror (_proc_,"Wrong particle statistics");

//   if (!WORM && BOSONS)
//   nrerror (_proc_,"BE statistics: worm algorithm only");

/*#ifndef HOSC_TEST 
   if (BOSONS && (IMTYPE < 0))  // need to define the reference axes for the area
   nrerror (_proc_,"Define the impurity type for the area estimator");
#endif */
}

void MCConfigInit(void)
{
   const char *_proc_=__func__;    // "MCConfigInit()";

#ifndef HOSC_TEST
   initLattice_config(MCCoords); 
   replInitial_config(MCCoords);

   if (NDIM != 3) nrerror(_proc_,"Only 3D for rotational coordinates");
#else
   for (int id=0;id<NDIM;id++)	
   for (int atom=0;atom<NumbAtoms;atom++)
   for (int it=0;it<NumbTimes;it++)
   MCCoords[id][atom*NumbTimes+it] = 0.0;	  
#endif

   cout<<"initial MCCoords "<<MCCoords[0][0]<<endl;

   for (int it=0;it<(NumbAtoms*NumbTimes);it++)
   {
       MCAngles[PHI][it] = 0.0;
       MCAngles[CTH][it] = 1.0;
//     toby
       MCAngles[CHI][it] = 0.0;

       double phi  = MCAngles[PHI][it];
      
       double cost = MCAngles[CTH][it];
       double sint = sqrt(1.0 - cost*cost);
  
       MCCosine[AXIS_X][it] = sint*cos(phi);
       MCCosine[AXIS_Y][it] = sint*sin(phi);
       MCCosine[AXIS_Z][it] = cost;
   }
}

void initLattice_config(double **pos)
// treat atoms and molecules separateley
// generate a cubic lattice if NumbAtoms = m^3, m - integer	
// nslices = Number of time slices	
{
   cout<<"in initLattice"<<endl;
   const char *_proc_ = __func__;    // "initLattice_config";

   int natoms = 0;                   // number of atoms 
   int nmolcs = 0;                   // number of molecules

   for (int type=0;type<NumbTypes;type++)
   if ((MCAtom[type].molecule == 1) || (MCAtom[type].molecule == 2)) nmolcs += MCAtom[type].numb; 
   else                       natoms += MCAtom[type].numb; 

// ----- INITIAL CONFIGURATION FOR ATOMS ----------------------

// box size per particle for atoms only: 
   double abox = BoxSize/pow((double)natoms,1.0/(double)NDIM); 
   double shift[NDIM];

   for (int id=0;id<NDIM;id++)
   shift[id] = 0.5*abox;  
    
   for (int type=0;type<NumbTypes;type++)   // count molecules only
   if (MCAtom[type].molecule == 0)
   {
      int offset = MCAtom[type].offset;         
      int maxnum = offset + MCAtom[type].numb*NumbTimes;         
      
      for (int atom=offset;atom<maxnum;atom+=NumbTimes)
      {
         for (int id=0;id<(NDIM-1);id++) 
         if (shift[id] > BoxSize) 
         {
            shift[id]    = 0.5*abox; 
            shift[id+1] += abox;
         }
 
         for (int id=0;id<NDIM;id++)   // set the center of the box at the origin
         pos[id][atom] = shift[id] - 0.5*BoxSize;

         shift[0] += abox;
      } 
   }

// ----- INITIAL CONFIGURATION FOR MOLECULES -------------------

// box size per particle for molecules only: 

   abox = BoxSize/pow((double)nmolcs,1.0/(double)NDIM); 

   cout<<"abox="<<abox<<" "<<BoxSize<<endl;
  
   for (int id=0;id<NDIM;id++)
   shift[id] = 0.5*abox;  
    
   for (int type=0;type<NumbTypes;type++)   // count molecules only
   if ((MCAtom[type].molecule == 1) || (MCAtom[type].molecule == 2) )
   {
      int offset = MCAtom[type].offset;         
      int maxnum = offset + MCAtom[type].numb*NumbTimes;         
      
      for (int atom=offset;atom<maxnum;atom+=NumbTimes)
      {  
        for (int id=0;id<(NDIM-1);id++) 
        if (shift[id] > BoxSize) {shift[id] = 0.5*abox; shift[id+1] += abox;}
 
        for (int id=0;id<NDIM;id++) // set the center of the box at the origin
        {
           pos[id][atom] = shift[id] - 0.5*BoxSize;

           if ((natoms == nmolcs))   // to avoid the overlap between particles
           pos[id][atom] += BoxSize;

           cout<<atom<<" "<<id<<" "<<pos[id][atom]<<endl;

        } 
 
        shift[0] += abox;
      }  // END loop over atoms
    }    // END loop over types
}

/*  no difference between atoms and molecules in the code below

void initLattice_config(double **pos, int nslices)
// generate a cubic lattice if NumbAtoms = m^3, m - integer	
// nslices = Number of time slices	
{
   const char *_proc_ = __func__;    // "initLattice_config";
 
// box size per particle: 
   double abox = BoxSize/pow((double)NumbAtoms,1.0/(double)NDIM); 
   double shift[NDIM];

   for (int id=0;id<NDIM;id++)
   shift[id] = 0.5*abox;  
    
   int max = NumbAtoms*nslices; 

   for (int atom=0;atom<max;atom+=nslices)
   {  
       for (int id=0;id<(NDIM-1);id++) 
       if (shift[id] > BoxSize) {shift[id] = 0.5*abox; shift[id+1] += abox;}
 
       for (int id=0;id<NDIM;id++)   // set the center of the box at the origin
       pos[id][atom] = shift[id] - 0.5*BoxSize;
 
       shift[0] += abox;
   }
}
*/

void replInitial_config(double **pos)
// replicate configurations for all time slices
{
   for (int id=0;id<NDIM;id++)	
   for (int atom=0;atom<NumbAtoms;atom++)
   for (int it=1;it<NumbTimes;it++)
   pos[id][atom*NumbTimes+it] = pos[id][atom*NumbTimes];	  
}	
