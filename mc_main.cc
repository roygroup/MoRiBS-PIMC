//
//         main()
//

#include <stdlib.h>
//#include <stdlib>
#include <sstream>
#include <iomanip>
#include <cmath>

#include "mc_confg.h"
#include "mc_setup.h"
#include "mc_input.h"
#include "mc_randg.h"
#include "mc_utils.h"
#include "mc_poten.h"
#include "mc_piqmc.h"
#include "mc_qworm.h"
#include "mc_estim.h"
#include "mc_const.h"

#include <stdio.h>
//#include <mpi.h>
#include <omp.h>
#include "rngstream.h"
#include "omprng.h"

void MCWormAverage(void);
void MCWormAverageReset(void);

double avergCount;   // # of calls of get_estim inside a block
double totalCount;   // sum avergCount   
         
void PIMCPass(int,int);

void MCGetAverage(void);

void MCResetBlockAverage(void);
void MCSaveBlockAverages(long int);

void MCSaveAcceptRatio(long int,long int,long int);

//--------------- BLOCK AVERAGE ------------

double _dbpot;       // potential energy differencies, block average  added by Hui Li
double _bpot;       // kinetic   energy, block average
double _bkin;       // potential energy, block average
double _bCv;        // heat capacity, block average
double _bCv_trans;  // translational heat capacity, block average
double _bCv_rot;   //  rotational heat capacity, block average

double _dpot_total;  // potential energy differences, global average  added by Hui Li
double _pot_total;  // kinetic   energy, global average
double _kin_total;  // potential energy, global average 

double _brot;       // rotational kin energy, block average
double _rot_total;  // rotational kin energy, global average
double _brotsq;     // rotational energy square, block average
double _rotsq_total; // rotational energy square, global average
double _Cv_total;    // heat capacity, global average
double _Cv_trans_total;    // translational heat capacity, global average
double _Cv_trans_1_total;    // translational heat capacity, global average
double _Cv_trans_2_total;    // translational heat capacity, global average
double _Cv_rot_total;    // rotational heat capacity, global average

fstream _feng;      // save accumulated energy

//---------------- ESTIMATORS ---------------

void SaveEnergy    (const char [],double,long int); // block average
void SaveSumEnergy (double,double);                 // accumulated average

void InitTotalAverage(void);
void DoneTotalAverage(void);

//---------------- INITIAL MCCOORDS -----------

extern "C" void initconf_(double *MCCooInit, double *MCAngInit, int *Ntotal, int *NBoson, int *PIndex,int *Rindex);

//---------------- PRINT PERMUTATION TABLE -----

extern "C" void prtper_(int *PIndex,int *NBoson,long int *blockCount);

//-------------------------------------------

int main(int argc, char *argv[])
{
//int numprocs, rank, namelen;
//char processor_name[MPI_MAX_PROCESSOR_NAME];


//MPI_Init(&argc, &argv);
//MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
//MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//MPI_Get_processor_name(processor_name, &namelen);
//MPI_Status mpistatus;

 
//-------------------------
   MPIsize = 1;            // this is for InitRandom()          
   MPIrank = MPI_MASTER;   // this is for InitRandom()          
//-----------------------------

   int restart = 0;
   IOReadParams(FINPUT,restart); // read input parameters 
// get the number of rotational steps treated by worker CPUs
// chunksize = NumbRotTimes / numprocs;

// pass numprocs to NProcs, which is a global variable
// NProcs = numprocs;

   MCInitParams();               // set system dependent parameters
   MCSetUnits(); 
   MCMemAlloc();

   MemAllocMCCounts();
   MemAllocQWCounts();

//head CPU comes in after memory declaration.
//if ( rank == MPI_MASTER ){// in the master cpu section

// only internal system of units after this point 
 
   MCInit();

   if (WORM)
   MCWormInit();

// clean the permutation sampling from the last run
   if (FileExist(FPERMU))
   _io_error("QMC ->",IO_ERR_FEXST,FPERMU);

   if (!restart) // new run, generate new status, rnd() streams and config files     
   {
// check the existence of the old files first
  
      if (FileExist(FSTATUS))    
     _io_error("QMC ->",IO_ERR_FEXST,FSTATUS);
       
      if (FileExist(FCONFIG)) 
     _io_error("QMC ->",IO_ERR_FEXST,FCONFIG);

      if (FileExist(FRANDOM)) 
     _io_error("QMC ->",IO_ERR_FEXST,FRANDOM);

//--------------------------------------------------------
//    generate tables - potentals, configurations etc, status

//    MCStartBlock = 0; 
//    SEED for head CPU
      SEED = 985456376;
      RandomInit(MPIrank,MPIsize);

      MCConfigInit();                // generate initial configurations

//    read in initial MCCoords and MCAngles
      if(InitMCCoords)
      {
         cout<<"read in MCCoords and MCAngles from xyz.init file"<<endl;
         int ntotal=NumbAtoms*NumbTimes;
         int iperm,nboson;
         nboson=0;
         if(BOSONS)
         {
            nboson=MCAtom[BSTYPE].numb;
         }
//       initconf_(MCCooInit,MCAngInit,&ntotal,&MCAtom[BSTYPE].numb,PIndex,RIndex);
         initconf_(MCCooInit,MCAngInit,&ntotal,&nboson,PIndex,RIndex);
         if(BOSONS)
         {
            cout<<"read PIndex: "<<" ";
            for(int id=0;id<MCAtom[BSTYPE].numb;id++)
            cout<<PIndex[id]<<" ";

            cout<<endl;

            cout<<"read RIndex: "<<" ";
            for(int id=0;id<MCAtom[BSTYPE].numb;id++)
            cout<<RIndex[id]<<" ";

            cout<<endl;
         }

         for (int it=0;it<NumbAtoms*NumbTimes;it++)
         {
            for(int id=0;id<NDIM;id++)
            {
               MCCoords[id][it]=MCCooInit[it*NDIM+id];
               MCAngles[id][it]=MCAngInit[it*NDIM+id];
            }

//          put into MCCosine
            double phi  = MCAngles[PHI][it];
            double cost = MCAngles[CTH][it];
            double sint = sqrt(1.0 - cost*cost);

            MCCosine[AXIS_X][it] = sint*cos(phi);
            MCCosine[AXIS_Y][it] = sint*sin(phi);
            MCCosine[AXIS_Z][it] = cost;


         }
      }

//    string fname = MCFileName + IO_EXT_XYZ;
//    IOxyz(IOWrite,fname.c_str());  // test output of initial config
      string fname = MCFileName;
      IOxyzAng(IOWrite,fname.c_str()); // test output of initial config

//--------------------------------------------------------

      StatusIO(IOWrite,FSTATUS);
      ConfigIO(IOWrite,FCONFIG);
      TablesIO(IOWrite,FTABLES);
      RandomIO(IOWrite,FRANDOM);

      if (WORM)
      QWormsIO(IOWrite,FQWORMS);
   }

   MCWormAverageReset();      // debug worm

// --- RESTART/START NEW RUN ----------------------------

   StatusIO(IORead,FSTATUS);  // load MCStartBlock
   ConfigIO(IORead,FCONFIG);  // load atoms/molecules positions
   TablesIO(IORead,FTABLES);  // load permutation tables
   RandomIO(IORead,FRANDOM);  // load rnd streams 

   if (WORM)
   QWormsIO(IORead,FQWORMS);

// BEGIN loop over blocks ------------------------

   InitPotentials();

   if (ROTATION)
   InitRotDensity();

   if((IREFLX == 1 || IREFLY == 1 || IREFLZ == 1) && MCAtom[IMTYPE].molecule != 2)
   {
      cout<<IREFLX<<" "<<IREFLY<<" "<<IREFLZ<<" "<<MCAtom[IMTYPE].molecule<<endl;
      nrerror("main","SYMMETRIZATION IS ONLY IMPLEMENTED FOR TOPS");
   }
// MPI BLOCK 1 in MASTER start
// send run information to all the worker CPUs
/*
   for(int iproc=1;iproc<numprocs;iproc++)
   {
      MPI_Send(&NumbTypes,1,MPI_INT,iproc,tagNumbTypes,MPI_COMM_WORLD);
      MPI_Send(&NumbAtoms,1,MPI_INT,iproc,tagNumbAtoms,MPI_COMM_WORLD);
      MPI_Send(&rhoprp[0],SizeRotDen,MPI_DOUBLE,iproc,tagrho,MPI_COMM_WORLD);
      MPI_Send(&erotpr[0],SizeRotDen,MPI_DOUBLE,iproc,tagero,MPI_COMM_WORLD);
      MPI_Send(&vtable[0],SizePotTab,MPI_DOUBLE,iproc,tagpot,MPI_COMM_WORLD);
      MPI_Send(&WORM,1,MPI_INT,iproc,tagWORM,MPI_COMM_WORLD);
      MPI_Send(&MCType[0],NumbAtoms,MPI_INT,iproc,tagMCType,MPI_COMM_WORLD);
      MPI_Send(&MCAtom[0].molecule,MAX_NUMBER_TYPES,MPI_INT,iproc,tagMCAtom_molecule,MPI_COMM_WORLD);
      MPI_Send(&RotRatio,1,MPI_INT,iproc,tagRotRatio,MPI_COMM_WORLD);
      MPI_Send(&Worm.type,1,MPI_INT,iproc,tagWormtype,MPI_COMM_WORLD);
      MPI_Send(&NumbRotTimes,1,MPI_INT,iproc,tagNumbRotTimes,MPI_COMM_WORLD);
      MPI_Send(&MCAtom[0].offset,MAX_NUMBER_TYPES,MPI_INT,iproc,tagMCAtom_offset,MPI_COMM_WORLD);
      MPI_Send(&IMTYPE,1,MPI_INT,iproc,tagIMTYPE,MPI_COMM_WORLD);
      MPI_Send(&MCAtom[0].numb,MAX_NUMBER_TYPES,MPI_INT,iproc,tagMCAtom_numb,MPI_COMM_WORLD);
      MPI_Send(&MCTau,1,MPI_DOUBLE,iproc,tagMCTau,MPI_COMM_WORLD);
      MPI_Send(&ROTATION,1,MPI_INT,iproc,tagROTATION,MPI_COMM_WORLD);
      MPI_Send(&MCStartBlock,1,MPI_LONG,iproc,tagMCStartBlock,MPI_COMM_WORLD);
      MPI_Send(&NumberOfMCBlocks,1,MPI_LONG,iproc,tagNumberOfMCBlocks,MPI_COMM_WORLD);
      MPI_Send(&NumberOfMCPasses,1,MPI_LONG,iproc,tagNumberOfMCPasses,MPI_COMM_WORLD);
//    MPI_Send(&TZMAT[1][1],NDIM,MPI_DOUBLE,iproc,tagTZMAT,MPI_COMM_WORLD);
      MPI_Send(&RotDenType,1,MPI_INT,iproc,tagRotDenType,MPI_COMM_WORLD);
      MPI_Send(&RotOdEvn,1,MPI_INT,iproc,tagRotOdEvn,MPI_COMM_WORLD);
      MPI_Send(&RotEoff,1,MPI_DOUBLE,iproc,tagRotEoff,MPI_COMM_WORLD);
      MPI_Send(&X_Rot,1,MPI_DOUBLE,iproc,tagARot,MPI_COMM_WORLD);
      MPI_Send(&Y_Rot,1,MPI_DOUBLE,iproc,tagBRot,MPI_COMM_WORLD);
      MPI_Send(&Z_Rot,1,MPI_DOUBLE,iproc,tagCRot,MPI_COMM_WORLD);
      MPI_Send(&MCRotTau,1,MPI_DOUBLE,iproc,tagMCRotTau,MPI_COMM_WORLD);
   }
// tagrunning = (tagNumberOfMCPasses + 2) % tagUpper;
// tagrunning = (tagTZMAT + 2) % tagUpper;
   tagrunning = (tagMCRotTau + 2) % tagUpper;

// the following loop passes the same tagrunning to all workers, in order to guarantee all CPUs have the same running tag
   for(int iproc=1;iproc<numprocs;iproc++)
   {
      MPI_Send(&tagrunning,1,MPI_INT,iproc,tagTAGRUN,MPI_COMM_WORLD);
   }
// MPI BLOCK 1 in MASTER done
*/

   InitMCEstims();
   InitTotalAverage();      // DUMP 
   double sumsCount = 0.0;  // counter for accum sums
   
   ResetMCCounts();
   ResetQWCounts();

/*
// try openmp for loop
   int chunksize,nthrds;
   #pragma omp parallel
   {
   int tid=omp_get_thread_num();
   nthrds=omp_get_num_threads();
   chunksize=NumbRotTimes/nthrds;
   cout<<"chunksize="<<chunksize<<endl;
   int itini=chunksize*tid;
   int itfnl=itini+chunksize;
// #pragma omp for
   for (int itrot=itini;itrot<itfnl-1;itrot++)
   {
      cout<<"itrot="<<itrot<<" tid="<<tid<<endl;
   }
   }

   for (int itrot=chunksize-1;itrot<NumbRotTimes;itrot=itrot+chunksize)
   {
      cout<<"itrot="<<itrot<<endl;
   }
   for (int itrot=nthrds*chunksize;itrot<NumbRotTimes;itrot++)
   {
      cout<<"itrot="<<itrot<<endl;
   }
*/

// get the chunksize and # of threads for the parallel rotation
   #pragma omp parallel
   {
      NThreads=omp_get_num_threads();
      chunksize=NumbRotTimes/NThreads;
   } // end omp

   cout<<"NThreads="<<NThreads<<" chunksize="<<chunksize<<endl;

   printf("OpenMP version: %d\n", _OPENMP);

   randomseed(); //set seed according to clock
// RngStream Rng[omp_get_num_procs()];     // initialize a parallel RNG named "Rng"

   long int blockCount = MCStartBlock;  
   while (blockCount<NumberOfMCBlocks+MCStartBlock) // START NEW BLOCK      
   {      
       blockCount++; 
       MCResetBlockAverage();
     
       long int passCount = 0;        // BEGIN NEW MC PASS
       long int passTotal = 0;        // total number of Moves = passCount*time 

       while (passCount++ < NumberOfMCPasses) 
       for (int time=0; time<NumbTimes; time++)
       {
           passTotal++;  

//         toby's printing
//         cout<<blockCount<<" "<<passCount<<" "<<time<<endl;
           for (int type=0;type<NumbTypes;type++)
           if (WORM && (type == Worm.type)) 
           {

              MCWormMove();

              if (!Worm.exists)
              {
/* reactive */
                 int rt = nrnd2(NumbTimes);   // select a time slice

                 if ((type == BSTYPE) || (type == FERMTYPE))
                 MCBisectionMoveExchange (type,rt);
                 else
                 MCBisectionMove (type,rt);
/* reactive */
#ifdef MOVECENTROIDTEST
                 if  (time == 0)
                 if ((type == BSTYPE) || (type == FERMTYPE))
                 MCMolecularMoveExchange(type);        
                 else
                 MCMolecularMove(type);        
#endif
              }
           } 
           else
           PIMCPass(type,time);

           if (blockCount>NumberOfEQBlocks)        // skip equilibration steps
           {
// evaluate averages 
              if (passTotal % MCSKIP_AVERG == 0)   // skip correlated configurations
              {
                 if (WORM)
                 {
                    if (!Worm.exists)
                    {
                       MCGetAverage();
//                     prtper_(PIndex,&MCAtom[BSTYPE].numb,&blockCount);

                       // print the instantaneous xyz and prl info for the closed path
                       if(PrintXYZprl)
                       {
                          stringstream bc;                // convert block # to string
                          bc.width(IO_BLOCKNUMB_WIDTH);
                          bc.fill('0');
                          bc<<blockCount;
                          string fname = MCFileName + bc.str();  // file name prefix including block #
                          IOxyzAng(IOWrite,fname.c_str());
                          PrintXYZprl = 0;
                       }

                    }
                    else  
                    MCWormAverage();  /// TEST ONLY
                 }   
                 else
                 {
                    MCGetAverage();
//                  omp_set_num_threads(1);
//                  MCGetAverage();
//                  omp_set_num_threads(NThreads);
                    // print the instantaneous xyz and prl info for the closed path
                    if(PrintXYZprl)
                    {
                       stringstream bc;                // convert block # to string
                       bc.width(IO_BLOCKNUMB_WIDTH);
                       bc.fill('0');
                       bc<<blockCount;
                       string fname = MCFileName + bc.str();  // file name prefix including block #
                       IOxyzAng(IOWrite,fname.c_str());
                       PrintXYZprl = 0;
                    }
                 }
              }

// DUMP, save global average. if avergCount = 0, then MCGetAverage is never called in this block.  All Saving steps are skipped.

              if (passTotal % MCSKIP_TOTAL == 0 && avergCount)  
              {
                 sumsCount += 1.0;                 
                 SaveSumEnergy (totalCount,sumsCount);
              }
           }  
          
           if (passTotal % MCSKIP_RATIO == 0)
           MCSaveAcceptRatio(passTotal,passCount,blockCount);
       }                               
       
// END  loop over MC passes (time slices) -----------------

       if (blockCount>NumberOfEQBlocks && avergCount)   // skip equilibration steps
       {
           MCSaveBlockAverages(blockCount);

//         save accumulated interatomic distribution
           SaveGraSum(MCFileName.c_str(),totalCount);

           if(IMPURITY && MCAtom[IMTYPE].molecule == 1)
           SaveDensities2D(MCFileName.c_str(),totalCount,MC_TOTAL);

           if(IMPURITY && MCAtom[IMTYPE].molecule == 2)
           {
              SaveDensities3D(MCFileName.c_str(),totalCount,MC_TOTAL); // this step takes lots of space.  temporarily turnned off
              SaveRho1D(MCFileName.c_str(),totalCount,MC_TOTAL);
//            SaveRhoThetaChi(MCFileName.c_str(),totalCount,MC_TOTAL); // we don't need 2d angular distribution comparison now
           }

         if (ROTATION)                  // DUMP  accumulated average
         SaveRCF(MCFileName.c_str(),totalCount,MC_TOTAL);
       }

//  CHECKPOINT: save status, rnd streams and configs ------

//     MCStartBlock = blockCount; 

       IOFileBackUp(FSTATUS); StatusIO(IOWrite,FSTATUS);
       IOFileBackUp(FCONFIG); ConfigIO(IOWrite,FCONFIG);
       IOFileBackUp(FTABLES); TablesIO(IOWrite,FTABLES);
       IOFileBackUp(FRANDOM); RandomIO(IOWrite,FRANDOM);      
      
      string fname = MCFileName + IO_EXT_XYZ;
      IOxyz(IOWrite,fname.c_str());  // test output of initial config

       if (WORM)
       {
          IOFileBackUp(FQWORMS);
          QWormsIO(IOWrite,FQWORMS);
       }
   }  
// END  loop over blocks ------------------------
// } // master cpu if block ends

// MPI_Finalize();

   DoneTotalAverage();  // DUMP

   if (WORM)
   MCWormDone();

   DoneMCEstims();
   DonePotentials();

   if (ROTATION)
   DoneRotDensity();

   MCMemFree();
   RandomFree();

   MFreeMCCounts();
   MFreeQWCounts();

   return 1; 
}

void PIMCPass(int type,int time)
{
   if (time == 0)
   MCMolecularMove(type);        

   MCBisectionMove(type,time);

   if ((type == IMTYPE) && ROTATION && MCAtom[type].molecule == 1)  // rotational degrees of freedom
   MCRotationsMove(type);

   if ((type == IMTYPE) && ROTATION && MCAtom[type].molecule == 2)  // non-linear rotor rotation added by Toby
   MCRotations3D(type);
}

void MCResetBlockAverage(void) 
{
   avergCount = 0.0;

   ResetMCEstims();

   ResetMCCounts();
   ResetQWCounts();

  _dbpot = 0.0;  //added by Hui Li
  _bpot = 0.0;
  _bkin = 0.0;

  _brot = 0.0;
  _brotsq = 0.0;
  _bCv    = 0.0;
  _bCv_trans = 0.0;
  _bCv_rot = 0.0;

   PrintXYZprl = 1;

   PrintYrfl = 1;
   PrintXrfl = 1;
   PrintZrfl = 1;

}

void MCGetAverage(void) 
{
   avergCount += 1.0;
   totalCount += 1.0;  

   double skin = GetKinEnergy();           // kin energy
   double spot = GetPotEnergy_Densities(); // pot energy and density distributions
   double srot = 0.0;
//   double dspot = GetPotEnergy_Diff(); // pot energy differencies added by Hui Li 

  _bkin += skin;  // block average for kin energy
  _bpot += spot;  // block average for pot energy
 // _dbpot += dspot;  // block average for pot energy differencies added by Hui Li
 
  _kin_total  += skin; // accumulated average 
  _pot_total  += spot;
//  _dpot_total  += dspot;  //added by Hui Li

// rotational degrees of freedom
/* reactive */
   if (ROTATION)
   {
//     double srot;

       if(MCAtom[IMTYPE].molecule == 1)
       srot = GetRotEnergy();           // kin energy
       else
       srot = GetRotE3D();
        
      _brot       += srot;
      _rot_total  += srot;

      _brotsq      += ErotSQ;
      _rotsq_total += ErotSQ;

       GetRCF(); 
   }

// accumulate terms for Cv
   double sCv;
   sCv = -0.5*(double)(NDIM*NumbAtoms)/(MCBeta*MCTau)
         -((double)(NDIM*NumbAtoms)*0.5/MCTau - skin - spot - srot)*((double)(NDIM*NumbAtoms)*0.5/MCTau - skin - spot - srot)
         + (2.0/MCBeta)*(0.5*(double)(NDIM*NumbAtoms)/MCTau - skin) + Erot_termSQ - ErotSQ;

   _bCv += sCv;
   _Cv_total += sCv;

// accumulate terms for translational Cv
   double sCv_trans;
   double sCv_trans_1=-((double)(NDIM*NumbAtoms)*0.5/MCTau - skin)*((double)(NDIM*NumbAtoms)*0.5/MCTau - skin);
   double sCv_trans_2=(2.0/MCBeta)*0.5*(double)(NDIM*NumbAtoms)/MCTau - skin;
   sCv_trans = -0.5*(double)(NDIM*NumbAtoms)/(MCBeta*MCTau)
         -((double)(NDIM*NumbAtoms)*0.5/MCTau - skin)*((double)(NDIM*NumbAtoms)*0.5/MCTau - skin)
         + (2.0/MCBeta)*(0.5*(double)(NDIM*NumbAtoms)/MCTau - skin);

   _bCv_trans += sCv_trans;
   _Cv_trans_total += sCv_trans;
   _Cv_trans_1_total += sCv_trans_1;
   _Cv_trans_2_total += sCv_trans_2;

// accumulate terms for rotational Cv
   double sCv_rot;
   sCv_rot = -srot*srot + Erot_termSQ - ErotSQ;

   _bCv_rot += sCv_rot;
   _Cv_rot_total += sCv_rot;

/* reactive */
// check whether there're bosons in the system

   if (BOSONS) 
   {
      GetExchangeLength();
      GetPermutation();
   }

#ifdef DEBUG_PIMC
   if (NDIM !=3) 
   nrerror("MCGetAverage()","Area estimators for 3D case only"); 
#endif
   if (BOSONS and MCAtom[IMTYPE].molecule == 1) 
   GetAreaEstimators();

   if (BOSONS)
   {
      int iframe = 0;
      GetAreaEstim3D(iframe);
   }

   if ((BOSONS and MCAtom[IMTYPE].molecule == 2) && ISPHER == 0) 
   {
      int iframe = 1;
      GetAreaEstim3D(iframe);
   }

//    reflect for MF molecule
      if(IREFLY == 1)
      {
         double rndrfl=rnd7();
         if(rndrfl < 0.5)
         {
            Reflect_MF_XZ();
         }
//       else
//       {
//          cout<<"Not in Reflect_MF_XZ"<<endl;
//       }
      }

      if(IREFLX == 1)
      {
         double rndrfl=rnd7();
         if(rndrfl < 0.5)
         {
            Reflect_MF_YZ();
         }
//       else
//       {
//          cout<<"Not in Reflect_MF_YZ"<<endl;
//       }
      }

      if(IREFLZ == 1)
      {
         double rndrfl=rnd7();
         if(rndrfl < 0.5)
         {
            Reflect_MF_XY();
         }
//       else
//       {
//          cout<<"Not in Reflect_MF_XY"<<endl;
//       }
      }

      if(IROTSYM == 1)
      {
         double rndrot=rnd7();

         if(rndrot < 0.5)
         RotSymConfig();
      }

}

void MCWormAverageReset(void)
{
}

void MCWormAverage(void)
{
}

void MCSaveBlockAverages(long int blocknumb) 
{
   const char *_proc_=__func__;    //   MCSaveBlockAverages()

   stringstream bc;                // convert block # to string
   bc.width(IO_BLOCKNUMB_WIDTH); 
   bc.fill('0');
   bc<<blocknumb;

   string fname = MCFileName + bc.str();  // file name prefix including block #

//-----------------------------------------------------------------
// densities

// toby temporarily by-pass the density save for the non-linear rotor
   if( IMPURITY && MCAtom[IMTYPE].molecule == 1)
   {
   SaveDensities1D    (fname.c_str(),avergCount);

   SaveDensities2D    (fname.c_str(),avergCount,MC_BLOCK);
   }

   if( IMPURITY && MCAtom[IMTYPE].molecule == 2)
   {
   SaveDensities1D    (fname.c_str(),avergCount);

   SaveRho1D          (fname.c_str(),avergCount,MC_BLOCK);

// SaveDensities3D    (fname.c_str(),avergCount,MC_BLOCK); // this step takes lots of space. temporarily turnned off

   SaveRhoThetaChi    (fname.c_str(),avergCount,MC_BLOCK);

// IOxyzAng(IOWrite,fname.c_str());
   }

   if (ROTATION) 
   SaveRCF            (fname.c_str(),avergCount,MC_BLOCK); 

   SaveEnergy         (MCFileName.c_str(),avergCount,blocknumb);

   if (BOSONS) 
   SaveExchangeLength (MCFileName.c_str(),avergCount,blocknumb);

   if (BOSONS && MCAtom[IMTYPE].molecule ==1)  
   SaveAreaEstimators (MCFileName.c_str(),avergCount,blocknumb);

   if (BOSONS)
   {
      int iframe = 0;
      SaveAreaEstim3D (MCFileName.c_str(),avergCount,blocknumb,iframe);
   }

   if ((BOSONS && MCAtom[IMTYPE].molecule ==2) && ISPHER == 0)
   {
      int iframe = 1;
      SaveAreaEstim3D (MCFileName.c_str(),avergCount,blocknumb,iframe);
   }

}

void SaveEnergy (const char fname [], double acount, long int blocknumb)
{
   const char *_proc_=__func__;    //  SaveEnergy()
 
   fstream fid;
   string fenergy;

   fenergy  = fname; 
   fenergy += IO_EXT_ENG; 

   fid.open(fenergy.c_str(),ios::app | ios::out);
   io_setout(fid);

   if (!fid.is_open())
  _io_error(_proc_,IO_ERR_FOPEN,fenergy.c_str());

   fid << setw(IO_WIDTH_BLOCK) << blocknumb  << BLANK;                 // block number
   fid << setw(IO_WIDTH) << _bkin*Units.energy/avergCount << BLANK;    // kinetic energy
   fid << setw(IO_WIDTH) << _bpot*Units.energy/avergCount << BLANK;    // potential anergy
   fid << setw(IO_WIDTH) <<(_bkin+_bpot)*Units.energy/avergCount << BLANK;   
   
// if (ROTATION)
   fid << setw(IO_WIDTH) << _brot*Units.energy/avergCount << BLANK;    // rot energy
   fid << setw(IO_WIDTH) << _brotsq*(Units.energy*Units.energy)/avergCount << BLANK;    // rot energy square
// fid << setw(IO_WIDTH) << _dbpot*Units.energy/avergCount << BLANK;    // potential differencies added by Hui Li 
   fid << setw(IO_WIDTH) <<(_bkin+_bpot+_brot)*Units.energy/avergCount << BLANK;  //total energy including rot energy 
   fid << setw(IO_WIDTH) << _bCv/avergCount << BLANK; // heat capacity
   fid << setw(IO_WIDTH) << _bCv_trans/avergCount << BLANK; // translational heat capacity
   fid << setw(IO_WIDTH) << _bCv_rot/avergCount << BLANK; // rotational heat capacity
   fid << endl;
   fid.close();
}

void SaveSumEnergy (double acount, double numb)  // global average
{
   const char *_proc_=__func__;    //  SaveSumEnergy()
 
  _feng << setw(IO_WIDTH_BLOCK) << numb << BLANK;    
  _feng << setw(IO_WIDTH) << _kin_total*Units.energy/acount << BLANK;    
  _feng << setw(IO_WIDTH) << _pot_total*Units.energy/acount << BLANK;    
  _feng << setw(IO_WIDTH) <<(_kin_total+_pot_total)*Units.energy/acount << BLANK;   

// if (ROTATION)
   {
  _feng << setw(IO_WIDTH) <<_rot_total*Units.energy/acount << BLANK;   
  _feng << setw(IO_WIDTH) <<_rotsq_total*(Units.energy*Units.energy)/acount << BLANK;   
   }
//_feng << setw(IO_WIDTH) << _dpot_total*Units.energy/acount << BLANK;   //added by Hui Li 
  _feng << setw(IO_WIDTH) <<(_kin_total+_pot_total+_rot_total)*Units.energy/acount << BLANK;  //total energy including rot  
// Cv
   double Cv = 0.5*(double)(NDIM*NumbAtoms*NumbTimes*Temperature)-(_kin_total+_pot_total+_rot_total)*Units.energy/acount;
   Cv = Cv*Cv + _Cv_total*Units.energy*Units.energy/acount;
   Cv = -Cv*MCBeta/Temperature;
  _feng << setw(IO_WIDTH) <<Cv << BLANK;

   double Cv_trans = 0.5*(double)(NDIM*NumbAtoms*NumbTimes*Temperature)-_kin_total*Units.energy/acount;
   Cv_trans = Cv_trans*Cv_trans+ _Cv_trans_total*Units.energy*Units.energy/acount;
   Cv_trans = -Cv_trans*MCBeta/Temperature;
  _feng << setw(IO_WIDTH) <<Cv_trans << BLANK;

// if (ROTATION)
   {
   double Cv_rot = -_rot_total*Units.energy/acount;
   Cv_rot = Cv_rot*Cv_rot+ _Cv_rot_total*Units.energy*Units.energy/acount;
   Cv_rot = -Cv_rot*MCBeta/Temperature;
  _feng << setw(IO_WIDTH) <<Cv_rot << BLANK;
   }

//_feng << setw(IO_WIDTH) << _Cv_trans_1_total*Units.energy*Units.energy/acount << BLANK;
//_feng << setw(IO_WIDTH) << _Cv_trans_2_total*Units.energy*Units.energy/acount << BLANK;

  _feng << endl;
}

void InitTotalAverage(void)  // DUMP
{
   const char *_proc_=__func__;    

   totalCount = 0.0;  // need to save in the status file 

  _kin_total = 0.0;   // need a function to reset all global average
  _pot_total = 0.0;
  _dpot_total = 0.0;  //added by Hui Li

  _rot_total = 0.0;
  _rotsq_total = 0.0;
  _Cv_total = 0.0;
  _Cv_trans_total = 0.0;
  _Cv_trans_1_total = 0.0;
  _Cv_trans_2_total = 0.0;
  _Cv_rot_total = 0.0;

// open files for output
// ENERGY

   string fenergy;

   fenergy  = MCFileName + IO_SUM; 
   fenergy += IO_EXT_ENG; 
 
   if (FileExist(fenergy.c_str()))   // backup the output of previous simulations 
   IOFileBackUp(fenergy.c_str());

  _feng.open(fenergy.c_str(), ios::out);
   io_setout(_feng);

   if (!_feng.is_open())
  _io_error(_proc_,IO_ERR_FOPEN,fenergy.c_str());
}

void DoneTotalAverage(void)
{
  _feng.close();
}

void MCSaveAcceptRatio(long int step,long int pass,long int block)
{
   int w = 8; 
 
   cout << "BLOCK:" << setw(w) << block << BLANK;
   cout << "PASS:"  << setw(w) << pass  << BLANK;
   cout << "STEP:"  << setw(w) << step  << BLANK;

   for (int type=0;type<NumbTypes;type++)
   if (WORM && (type == Worm.type))
   {
      double ratio_open    = QWAccep[0][QW_OPEN] /QWTotal[0][QW_OPEN];
      double ratio_close   = QWAccep[0][QW_CLOSE]/QWTotal[0][QW_CLOSE];

      double ratio_advance = QWAccep[0][QW_ADVANCE]/QWTotal[0][QW_ADVANCE];
      double ratio_recede  = QWAccep[0][QW_RECEDE] /QWTotal[0][QW_RECEDE];

      double ratio_swap    = QWAccep[0][QW_SWAP]/QWTotal[0][QW_SWAP];

      cout<<setw(w)<<"open/close";
      cout<<" [ "; 
      cout<<QWTotal[0][QW_OPEN]/countQW<<"-"<<QWTotal[0][QW_CLOSE]/countQW; 
      cout<<" ] "; 
      cout<<BLANK;                     
 
      cout<<setw(w)<<ratio_open   << BLANK; // accept ratio for "open"   move
      cout<<setw(w)<<ratio_close  << BLANK; // accept ratio for "close"  move

      cout<<setw(w)<<"advance/recede";
      cout<<" [ "; 
      cout<<QWTotal[0][QW_ADVANCE]/countQW<<"-"<<QWTotal[0][QW_RECEDE]/countQW; 
      cout<<" ] "; 
      cout<<BLANK;                     
 
      cout<<setw(w)<<ratio_advance<< BLANK;  // accept ratio for "advance" move 
      cout<<setw(w)<<ratio_recede << BLANK;  // accept ratio for "recede"  move

      cout<<setw(w)<<"swap";
      cout<<" [ "; 
      cout<<QWTotal[0][QW_SWAP]/countQW; 
      cout<<" ] "; 
      cout<<BLANK;                     
      cout<<setw(w)<<ratio_swap << BLANK;  // accept ratio for "recede"  move
   }
   else
   {
      double ratio_molec = MCAccep[type][MCMOLEC]/MCTotal[type][MCMOLEC];
      double ratio_multi = MCAccep[type][MCMULTI]/MCTotal[type][MCMULTI];
 
      cout<<setw(w)<<MCAtom[type].type<<BLANK; // atom type

      cout<<setw(w)<<ratio_molec<<BLANK;       // accept ratio for "molecular" move
      cout<<setw(w)<<ratio_multi<<BLANK;       // accept ratio for multilevel move 
   }

   if (ROTATION)
   {
      double ratio_rotat = MCAccep[IMTYPE][MCROTAT]/MCTotal[IMTYPE][MCROTAT];
      cout<<BLANK<< "Rot: ";
      cout<<setw(w)<<ratio_rotat<<BLANK;       // accept ratio for "molecular" move
   }

   cout << endl;

}
