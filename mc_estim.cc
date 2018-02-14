// estimators for MC

#include <cmath>
#include <iomanip>

//#include <gsl/gsl_sf_legendre.h>

#include "mc_setup.h"
#include "mc_poten.h"
#include "mc_utils.h"
#include "mc_estim.h"
#include "mc_input.h"
#include "mc_qworm.h"
#include "mc_confg.h"
#include "mc_const.h"

#include <omp.h>

// -----  density distributions ----------------------------

int NUMB_DENS1D;                // number of 1D density distributions
int NUMB_DENS2D;                // number of 2D density distributions
int NUMB_DENS3D;                // number of 3D density distributions added by Toby

const int MC_BINSR = 300;       // number of bins for radius original value=500
const int MC_BINST = 50;       // number of bins for theta original value=100
const int MC_BINSC = 100;       // number of bins for chi, added by Toby original value=200

const double MIN_RADIUS = 0.0;  // [AA], min radius for radial distribution functions 
const double MAX_RADIUS = 15.0; // [AA], max radius for radial distribution functions original value=25

double _delta_radius;           // \delta r for density distributions    
double _delta_theta;            // \delta\theta density distributions 
double _delta_chi;              // \delta\chi   density distributions

double _max_radius;             // max radius for density distributions
double _min_radius;             // min radius for density distributions 

double **   _gr1D;              // 1D radial distribution functions
double **   _gr1D_sum;          // 1D radial distribution functions accumulated
double ***  _gr2D;              // 2D distribution functions
double ***  _gr2D_sum;          // 2D distribution functions accumulated

double **   _gr3D;              // 3D distribution functions, added by Toby
double **   _gr3D_sum;          // 3D distribution functions accumulated, added by Toby

double *    _relthe_sum;        // Binning of relative euler angle theta between two rotational imaginary slices
double *    _relphi_sum;        // Binning of relative euler angle phi between two rotational imaginary slices
double *    _relchi_sum;        // Binning of relative euler angle chi between two rotational imaginary slices

double *    _drgid;             // radial  grid for distribution functions
double *    _trgid;             // angular grid for distribution functions

void densities_malloc(void);
void densities_init  (void);
void densities_mfree (void);
void densities_reset (int);  //revised by Hui Li

void bin_1Ddensity(double,int);
void bin_2Ddensity(double,double,int);
void bin_3Ddensity(double,double,double,int);

// -----  Rotational correlation functions ----

double ** _rcf;       
double ** _rcf_sum;

void rcf_malloc(void);
void rcf_init  (void);
void rcf_mfree (void);
void rcf_reset (int);

// SUPER (area estimator)

double _areas[2];
double _area2[2];

double _inert[2];

const int PERP = 0;
const int PARL = 1;

// SUPER (area estimator) for non-linear dopant
double _areas3DMFF[6];   // block accumulated area tensor estimator in dopant-fixed frame
double _inert3DMFF[9];  // block accumulated classical moment of inertia in dopant-fixed frame

double _areas3DSFF[6];   // block accumulated area tensor estimator in space-fixed frame
double _inert3DSFF[9];  // block accumulated classical moment of inertia in space-fixed frame

void super_reset (void);

// integer flag
int PrintXYZprl; // for printing instantaneous xyz coordinates and permutation table for each block

// rotational energy square estimator
double ErotSQ; // asymmetric top rotational energy square estimator
double Erot_termSQ; // sum of square of each bead's rotational energy estimator


//----------------------------------------------------

int     * _pflags;
double  * _ploops;

extern "C" void rotden_(double *Eulan1,double *Eulan2,double *Eulrel,double *rho,double *erot,double *esq,double *rhoprp, double *erotpr, double *erotsq, int *istop); // external fortran by toby

extern "C" void vcord_(double *Eulang, double *RCOM, double *Rpt, double *vtable, int *Rgrd,int *THgrd, int *CHgrd, double *Rvmax, double *Rvmin, double *Rvstep, double *vpot3d, double *radret, double *theret, double *chiret, double *hatx, double *haty, double *hatz, int *ivcord);

extern "C" void rsrot_(double *Eulan1,double *Eulan2,double *X_Rot,double *Y_Rot,double *Z_Rot,double *tau,int *iodevn,double *eoff,double *rho,double *erot);

extern "C" void rsline_(double *X_Rot,double *p0,double *tau,double *rho,double *erot);

extern "C" void vspher_(double *r,double *v);

extern "C" void caleng_(double *com_1, double *com_2, double *E_2H2O, double *Eulang_1, double *Eulang_2);

void InitMCEstims(void)
{
   densities_init();   
   densities_malloc();
   densities_reset (MC_TOTAL);  //revised by Hui 

  _pflags = new int    [MCAtom[BSTYPE].numb];
  _ploops = new double [MCAtom[BSTYPE].numb];

   if (ROTATION)
   {
      rcf_init();
      rcf_malloc();
      rcf_reset(MC_TOTAL);
   }

// SAVE RESULTS OF PREVIOUS SIMULATIONS

   string fname;

// energy

   fname  = MCFileName.c_str(); 
   fname += IO_EXT_ENG; 

   if (FileExist(fname.c_str()))
   { 
      IOFileBackUp(fname.c_str());
      IOFileDelete(fname.c_str());
   }

// area estimators

   fname  = MCFileName.c_str(); 
   fname += IO_EXT_SUP; 

   if (FileExist(fname.c_str()))
   { 
      IOFileBackUp(fname.c_str());
      IOFileDelete(fname.c_str());
   }

// 3d area estimators in dopant-fixed frame

   fname  = MCFileName.c_str();
   fname += IO_EXT_MFFSUP3D;

   if (FileExist(fname.c_str()))
   {
      IOFileBackUp(fname.c_str());
      IOFileDelete(fname.c_str());
   }

// 3d area estimators in space-fixed frame

   fname  = MCFileName.c_str();
   fname += IO_EXT_SFFSUP3D;

   if (FileExist(fname.c_str()))
   {
      IOFileBackUp(fname.c_str());
      IOFileDelete(fname.c_str());
   }

// exchange length

   fname  = MCFileName.c_str(); 
   fname += IO_EXT_PRL; 

   if (FileExist(fname.c_str()))
   { 
      IOFileBackUp(fname.c_str());
      IOFileDelete(fname.c_str());
   }
}

void DoneMCEstims(void)
{
   densities_mfree();

   delete [] _pflags;
   delete [] _ploops;

   if (ROTATION)
   rcf_mfree();
}

void ResetMCEstims(void)  // reset BLOCK averages
{
   densities_reset(MC_BLOCK);  //revised by Hui

   if (BOSONS)
   for (int atom=0;atom<MCAtom[BSTYPE].numb;atom++)
   {
      _pflags [atom] = 0;
      _ploops [atom] = 0.0;
   }

   if (ROTATION)
   rcf_reset(MC_BLOCK);

   super_reset();
}

// -----  BEGIN densities ------------------------------

void densities_init(void)
//
// atom-atom      (the same type) and 
// atom-molecule  distribution functions  
//
// LIMITATION:    the density type corresponds to the atom type
//                no cross densities  
{
   const char *_proc_=__func__;    // "densities_init()";

   if ((NUMB_ATOMTYPES > 1) || (NUMB_MOLCTYPES > 1))
   nrerror(_proc_,"No more then one atom/molecule type: densities and potential energy");   

// if ((NUMB_MOLCTYPES>0) && (NUMB_MOLCS!=1))
// nrerror(_proc_,"Only one molecular impurity");

   NUMB_DENS1D = NUMB_ATOMTYPES;  //# atom-atom densities [no cross-distributions]

   if (IMPURITY && (MCAtom[IMTYPE].molecule == 1))                          
   NUMB_DENS2D = NUMB_ATOMTYPES;  //# molecule-atoms distributions 

// Toby adds 3d distribution
   if (IMPURITY && (MCAtom[IMTYPE].molecule == 2))                          
   NUMB_DENS3D = NUMB_ATOMTYPES+NUMB_MOLCTYPES;  //# non-linear molecule-atoms distributions 

  _max_radius   = MAX_RADIUS/Units.length;  // max radius for radial distributions;
  _min_radius   = MIN_RADIUS/Units.length;  // min radius for radial distributions;

  _delta_radius = (_max_radius - _min_radius)/(double)MC_BINSR;
  _delta_theta  =  M_PI/(double)(MC_BINST - 1);
// Toby adds chi bin
  _delta_chi    =  2.0*M_PI/(double)(MC_BINSC - 1);
}

void densities_malloc(void)
// memory allocation for densities
{
  _gr1D  = doubleMatrix(NUMB_DENS1D,MC_BINSR);
  _gr1D_sum  = doubleMatrix(NUMB_DENS1D,MC_BINSR);

   if (IMPURITY && (MCAtom[IMTYPE].molecule == 1))
   {
      _gr2D  = new double ** [NUMB_DENS2D];
      _gr2D_sum  = new double ** [NUMB_DENS2D];  //added by Hui Li

      for (int id=0;id<NUMB_DENS2D;id++) 
      {
      _gr2D[id] = doubleMatrix(MC_BINSR,MC_BINST);
      _gr2D_sum[id] = doubleMatrix(MC_BINSR,MC_BINST); //added by Hui Li
      }
   }

// the following if block added by Toby
  if (IMPURITY && (MCAtom[IMTYPE].molecule == 2))
  {
     _gr3D  = doubleMatrix(NUMB_DENS3D,MC_BINSR*MC_BINST*MC_BINSC);
     _gr3D_sum  = doubleMatrix(NUMB_DENS3D,MC_BINSR*MC_BINST*MC_BINSC);

     _relthe_sum = new double [MC_BINST];
     _relphi_sum = new double [MC_BINSC];
     _relchi_sum = new double [MC_BINSC];

  }

}

void densities_mfree(void)
{
   free_doubleMatrix(_gr1D);
   free_doubleMatrix(_gr1D_sum);
   
   if (IMPURITY && MCAtom[IMTYPE].molecule == 1)
   {
      for (int id=0;id<NUMB_DENS2D;id++) 
      {
      free_doubleMatrix(_gr2D[id]);
      free_doubleMatrix(_gr2D_sum[id]);    //added by Hui Li
      }
      delete [] _gr2D;
      delete [] _gr2D_sum;   //added by Hui Li
    }

   if (IMPURITY && MCAtom[IMTYPE].molecule == 2)
   {
      free_doubleMatrix(_gr3D);
      free_doubleMatrix(_gr3D_sum);
   }

}

void densities_reset(int mode)     //revised by Hui Li
{
#ifdef DEBUG_PIMC
   const char *_proc_= __func__; //  rcf_reset() 

   if ((mode != MC_BLOCK) && (mode != MC_TOTAL))
   nrerror(_proc_,"Unknow mode");
#endif

    for (int id=0;id<NUMB_DENS1D;id++) 
    for (int ir=0;ir<MC_BINSR;ir++) 
    {
       if(mode == MC_BLOCK)
       _gr1D[id][ir] = 0.0;

       if(mode == MC_TOTAL)
       _gr1D_sum[id][ir] = 0.0;
    }

    if (IMPURITY && (MCAtom[IMTYPE].molecule == 1))
    for (int id=0;id<NUMB_DENS2D;id++) 
    for (int ir=0;ir<MC_BINSR;ir++) 
    for (int it=0;it<MC_BINST;it++) 
    {
       if(mode == MC_BLOCK)
       _gr2D[id][ir][it] = 0.0;     // block average

       if(mode == MC_TOTAL)        
       _gr2D_sum[id][ir][it] = 0.0;     // accumulated average
    }
    if (IMPURITY && (MCAtom[IMTYPE].molecule == 2))
    {
    for (int id=0;id<NUMB_DENS3D;id++)
    for (int ir=0;ir<MC_BINSR;ir++)
    for (int it=0;it<MC_BINST;it++)
    for (int ic=0;ic<MC_BINSC;ic++)
    {

       int ijk = (ir*MC_BINST + it)*MC_BINSC + ic;

       if(mode == MC_BLOCK)
       _gr3D[id][ijk] = 0.0;     // block average

       if(mode == MC_TOTAL)
       _gr3D_sum[id][ijk] = 0.0;     // accumulated average

    }

    if(mode == MC_TOTAL)
    {

       for (int it=0;it<MC_BINST;it++)
       _relthe_sum[it]=0.0;

       for (int ic=0;ic<MC_BINSC;ic++)
       {
          _relphi_sum[ic]=0.0;
          _relchi_sum[ic]=0.0;
       }

    }

    }
}

//------- RCF -------------------

void rcf_init(void)
{
}

void rcf_malloc(void)
// memory allocation for rotational correlation functions
{
  _rcf     = doubleMatrix(NUMB_RCF,NumbRotTimes);
  _rcf_sum = doubleMatrix(NUMB_RCF,NumbRotTimes);
}

void rcf_mfree(void)
{
   free_doubleMatrix(_rcf);
   free_doubleMatrix(_rcf_sum);
}

void rcf_reset(int mode)
{
#ifdef DEBUG_PIMC
   const char *_proc_= __func__; //  rcf_reset() 

   if ((mode != MC_BLOCK) && (mode != MC_TOTAL))
   nrerror(_proc_,"Unknow mode"); 
#endif 

   for (int ip=0;ip<NUMB_RCF; ip++) 
   for (int it=0;it<NumbRotTimes;it++)
   { 
      _rcf    [ip][it] = 0.0;  // block average

       if (mode == MC_TOTAL)   // total average
      _rcf_sum[ip][it] = 0.0;
    }
}
/*
//added by Hui Li
double GetPotEnergy_Diff(void)
{
   const char *_proc_=__func__; //  GetPotEnergy_Diff()  

#ifdef DEBUG_WORM
   if (Worm.exists)
   nrerror(_proc_," Only for Z-configurations");
#endif

   double dr[NDIM];
   double dspot = 0.0;

   for (int atom0=0;atom0<(NumbAtoms-1);atom0++)      
   for (int atom1=(atom0+1);atom1<NumbAtoms;atom1++)
   {
      int type0   = MCType[atom0];
      int type1   = MCType[atom1];

      int offset0 = NumbTimes*atom0;
      int offset1 = NumbTimes*atom1;

      for (int it=0;it<NumbTimes;it++) 	    
      {  
         int t0 = offset0 + it;
         int t1 = offset1 + it;

         double dr2 = 0.0;  		 
         for (int id=0;id<NDIM;id++)
         {
            dr[id]  = (MCCoords[id][t0] - MCCoords[id][t1]);

            if (MINIMAGE)
            dr[id] -= (BoxSize[id]*rint(dr[id]/BoxSize[id]));

            dr2    += (dr[id]*dr[id]);
         }
   	 
//#ifdef _CUTOFF_	     
//       if (dr2<dljcutoff2)
//#endif
         double r = sqrt(dr2);

//----------- [ATOM - MOLECULE] ----------------------

         if (MCAtom[type0].molecule||MCAtom[type1].molecule)  // 2D interaction 
         {
         //  type 1 is a molecule 
     
             int sgn   = 1;             // set to -1 to correct the orientaion of dr

             int tm    = offset1 + it/RotRatio;
//           int tm    = offset1 + floor((double)it/(double)RotRatio);

             int typep = type1;         // define the type of the potential
             int typed = type0;         // define the type of the density

         //  type 0 is a molecule ?   
             if (MCAtom[type0].molecule)  // does not work for two molecules
             {
                sgn   = -1;   
                tm    = offset0 + it/RotRatio;
                typep = type0; 
                typed = type1; 
             }

             double cost = 0.0;
             for (int id=0;id<NDIM;id++)    // n*dr = r*cos(theta) 
             cost += (MCCosine[id][tm]*dr[id]);   	 
	 
             cost /= r;                     // cos(theta)
             cost *= sgn;                   // correct the orientation 

             dspot += DLPot2D(r,cost,typep);  // potential differencies  
         }
      }  // LOOP OVER TIME SLICES
   }     // LOOP OVER ATOM PAIRS

   return (dspot/(double)NumbTimes);
}
*/



double GetPotEnergy_Densities(void)
// should be compatible with PotEnergy() from mc_piqmc.cc
{
   const char *_proc_=__func__; //  GetPotEnergy_Densities()  

#ifdef DEBUG_WORM
   if (Worm.exists)
   nrerror(_proc_," Only for Z-configurations");
#endif

// double dr[NDIM];
   double spot = 0.0;

   for (int atom0=0;atom0<(NumbAtoms-1);atom0++)      
   for (int atom1=(atom0+1);atom1<NumbAtoms;atom1++)
   {
      int type0   = MCType[atom0];
      int type1   = MCType[atom1];

      int offset0 = NumbTimes*atom0;
      int offset1 = NumbTimes*atom1;

      double spot_pair=0.0;

      #pragma omp parallel for reduction(+: spot_pair)
      for (int it=0;it<NumbTimes;it++) 	    
      {  
         double dr[NDIM];
         int t0 = offset0 + it;
         int t1 = offset1 + it;

         double dr2 = 0.0;  		 
         for (int id=0;id<NDIM;id++)
         {
            dr[id]  = (MCCoords[id][t0] - MCCoords[id][t1]);

            if (MINIMAGE)
            dr[id] -= (BoxSize[id]*rint(dr[id]/BoxSize[id]));

            dr2    += (dr[id]*dr[id]);
         }
   	 
//#ifdef _CUTOFF_	     
//       if (dr2<dljcutoff2)
//#endif
         double r = sqrt(dr2);

//----------- [ATOM - MOLECULE] ----------------------

         if ((MCAtom[type0].molecule == 1)||(MCAtom[type1].molecule == 1))  // 2D interaction 
         {
         //  type 1 is a molecule 
     
             int sgn   = 1;             // set to -1 to correct the orientaion of dr

             int tm    = offset1 + it/RotRatio;
//           int tm    = offset1 + floor((double)it/(double)RotRatio);

             int typep = type1;         // define the type of the potential
             int typed = type0;         // define the type of the density

         //  type 0 is a molecule ?   
             if (MCAtom[type0].molecule == 1)  // does not work for two molecules
             {
                sgn   = -1;   
                tm    = offset0 + it/RotRatio;
                typep = type0; 
                typed = type1; 
             }

             double cost = 0.0;
             for (int id=0;id<NDIM;id++)    // n*dr = r*cos(theta) 
             cost += (MCCosine[id][tm]*dr[id]);   	 
	 
             cost /= r;                     // cos(theta)
             cost *= sgn;                   // correct the orientation 

             bin_2Ddensity (r,cost,typed);  // densities 
             spot_pair += LPot2D(r,cost,typep);  // potential energy 
//           cout<<"it="<<it<<" spot_pair="<<spot_pair<<" LPot2D="<<LPot2D(r,cost,typep)<<endl;
        }
//----------- [ATOM - NON-LINEAR MOLECULE] ----------------------
        else if (((MCAtom[type0].molecule == 2)||(MCAtom[type1].molecule == 2)) && (MCAtom[type0].molecule != MCAtom[type1].molecule) ) // 3D interaction, no density is calculated now
        {

            int tm;
            int typed;

            double RCOM[3];
            double Rpt[3];
            double Eulang[3];
            double vpot3d;
            double radret;
            double theret;
            double chiret;
            double hatx[3];
            double haty[3];
            double hatz[3];
            int    ivcord = 0;
            if(MCAtom[type0].molecule == 2)
            {
//             determine type of atoms for bin_3Ddensity
               typed = type1;
               tm  = offset0 + it/RotRatio;
               for (int id=0;id<NDIM;id++)
               {
                  RCOM[id] = MCCoords[id][t0];
                  Rpt[id]  = MCCoords[id][t1];
               }
            }
            else
            {
//             determine type of atoms for bin_3Ddensity
               typed = type0;
               tm  = offset1 + it/RotRatio;
               for (int id=0;id<NDIM;id++)
               {
                  Rpt[id]  = MCCoords[id][t0];
                  RCOM[id] = MCCoords[id][t1];
               }
            }
            Eulang[PHI]=MCAngles[PHI][tm];
            Eulang[CTH]=acos(MCAngles[CTH][tm]);
            Eulang[CHI]=MCAngles[CHI][tm];

            if( ISPHER == 0)
            {
               vcord_(Eulang,RCOM,Rpt,vtable,&Rgrd,&THgrd,&CHgrd,&Rvmax,&Rvmin,&Rvstep,&vpot3d,&radret,&theret,&chiret,hatx,haty,hatz,&ivcord);
            }
            else if( ISPHER == 1)
            {
               radret = r;
               vspher_(&radret,&vpot3d);
               theret = 0.0;
               chiret = 0.0;
            }

            bin_3Ddensity (radret,theret,chiret,typed);  // accumulate density

            spot_pair += vpot3d;
        }
//---------[NON-LINEAR - NON-LINEAR from GG]
         else if ( ((MCAtom[type0].molecule == 2) && (MCAtom[type1].molecule == 2)) && (MCAtom[IMTYPE].numb > 1) )
         {
          //   if ( (MCType[atom0] == IMTYPE) && (MCType[atom1] == IMTYPE) )
         //   {
         //     if ( (MCAtom[type0].molecule == 2) && (MCAtom[type1].molecule == 2) )
         //   {
        //     cout<<"Le if de GG: MCAtom[type0].molecule MCAtom[type1].molecule"<<MCAtom[type0].molecule<<" "<<MCAtom[type1].molecule<<endl;
             double com_1[3];
             double com_2[3];
             double Eulang_1[3];
             double Eulang_2[3];
             double E_2H2O;
             for (int id=0;id<NDIM;id++)
             {
                  com_1[id] = MCCoords[id][t0];
                  com_2[id] = MCCoords[id][t1];
             }
             int tm0=offset0 + it/RotRatio;
             int tm1=offset1 + it/RotRatio;
             Eulang_1[PHI]=MCAngles[PHI][tm0];
             Eulang_1[CTH]=acos(MCAngles[CTH][tm0]);
             Eulang_1[CHI]=MCAngles[CHI][tm0];
             Eulang_2[PHI]=MCAngles[PHI][tm1];
             Eulang_2[CTH]=acos(MCAngles[CTH][tm1]);
             Eulang_2[CHI]=MCAngles[CHI][tm1];
             caleng_(com_1, com_2, &E_2H2O,
                        Eulang_1, Eulang_2);
             spot_pair += E_2H2O;
          //  }  //
          //  }
         }
//------------- [ATOM - ATOM] ------------------------------- 
         else                // only one atom type    
         if ((type0 == type1) && MCAtom[type0].molecule == 0) // no "cross" densities 
         {
            bin_1Ddensity (r,type1);    // densities 
            spot_pair += SPot1D(r,type1);    // potential energy
         }
      }  // LOOP OVER TIME SLICES
      spot += spot_pair;
   }     // LOOP OVER ATOM PAIRS

// cout<<"in GetPotDensity"<<" _gr1D[0][80]="<<_gr1D[0][80]<<" _gr1D_sum[0][80]="<<_gr1D_sum[0][80]<<endl;
// cout<<"spot="<<spot<<endl;
   return (spot/(double)NumbTimes);
}

double GetPotEnergy(void)
// should be compatible with PotEnergy() from mc_piqmc.cc
{
   const char *_proc_=__func__; //  GetPotEnergy_Densities()  

#ifdef DEBUG_WORM
   if (Worm.exists)
   nrerror(_proc_," Only for Z-configurations");
#endif

// double dr[NDIM];
   double spot = 0.0;

   for (int atom0=0;atom0<(NumbAtoms-1);atom0++)      
   for (int atom1=(atom0+1);atom1<NumbAtoms;atom1++)
   {
      int type0   = MCType[atom0];
      int type1   = MCType[atom1];

      int offset0 = NumbTimes*atom0;
      int offset1 = NumbTimes*atom1;

      double spot_pair=0.0;

      #pragma omp parallel for reduction(+: spot_pair)
      for (int it=0;it<NumbTimes;it++) 	    
      {  
         double dr[NDIM];
         int t0 = offset0 + it;
         int t1 = offset1 + it;

         double dr2 = 0.0;  		 
         for (int id=0;id<NDIM;id++)
         {
            dr[id]  = (MCCoords[id][t0] - MCCoords[id][t1]);

            if (MINIMAGE)
            dr[id] -= (BoxSize[id]*rint(dr[id]/BoxSize[id]));

            dr2    += (dr[id]*dr[id]);
         }
   	 
//#ifdef _CUTOFF_	     
//       if (dr2<dljcutoff2)
//#endif
         double r = sqrt(dr2);

//----------- [ATOM - MOLECULE] ----------------------

         if ((MCAtom[type0].molecule == 1)||(MCAtom[type1].molecule == 1))  // 2D interaction 
         {
         //  type 1 is a molecule 
     
             int sgn   = 1;             // set to -1 to correct the orientaion of dr

             int tm    = offset1 + it/RotRatio;
//           int tm    = offset1 + floor((double)it/(double)RotRatio);

             int typep = type1;         // define the type of the potential
             int typed = type0;         // define the type of the density

         //  type 0 is a molecule ?   
             if (MCAtom[type0].molecule == 1)  // does not work for two molecules
             {
                sgn   = -1;   
                tm    = offset0 + it/RotRatio;
                typep = type0; 
                typed = type1; 
             }

             double cost = 0.0;
             for (int id=0;id<NDIM;id++)    // n*dr = r*cos(theta) 
             cost += (MCCosine[id][tm]*dr[id]);   	 
	 
             cost /= r;                     // cos(theta)
             cost *= sgn;                   // correct the orientation 

//           bin_2Ddensity (r,cost,typed);  // densities 
             spot_pair += LPot2D(r,cost,typep);  // potential energy 
        }
//----------- [ATOM - NON-LINEAR MOLECULE] ----------------------
        else if (((MCAtom[type0].molecule == 2)||(MCAtom[type1].molecule == 2)) && (MCAtom[type0].molecule != MCAtom[type1].molecule) ) // 3D interaction, no density is calculated now
        {

            int tm;
            int typed;

            double RCOM[3];
            double Rpt[3];
            double Eulang[3];
            double vpot3d;
            double radret;
            double theret;
            double chiret;
            double hatx[3];
            double haty[3];
            double hatz[3];
            int    ivcord = 0;
            if(MCAtom[type0].molecule == 2)
            {
//             determine type of atoms for bin_3Ddensity
               typed = type1;
               tm  = offset0 + it/RotRatio;
               for (int id=0;id<NDIM;id++)
               {
                  RCOM[id] = MCCoords[id][t0];
                  Rpt[id]  = MCCoords[id][t1];
               }
            }
            else
            {
//             determine type of atoms for bin_3Ddensity
               typed = type0;
               tm  = offset1 + it/RotRatio;
               for (int id=0;id<NDIM;id++)
               {
                  Rpt[id]  = MCCoords[id][t0];
                  RCOM[id] = MCCoords[id][t1];
               }
            }
            Eulang[PHI]=MCAngles[PHI][tm];
            Eulang[CTH]=acos(MCAngles[CTH][tm]);
            Eulang[CHI]=MCAngles[CHI][tm];

            if( ISPHER == 0)
            {
               vcord_(Eulang,RCOM,Rpt,vtable,&Rgrd,&THgrd,&CHgrd,&Rvmax,&Rvmin,&Rvstep,&vpot3d,&radret,&theret,&chiret,hatx,haty,hatz,&ivcord);
            }
            else if( ISPHER == 1)
            {
               radret = r;
               vspher_(&radret,&vpot3d);
               theret = 0.0;
               chiret = 0.0;
            }

//          bin_3Ddensity (radret,theret,chiret,typed);  // accumulate density

            spot_pair += vpot3d;
        }
//---------[NON-LINEAR - NON-LINEAR from GG]
         else if ( ((MCAtom[type0].molecule == 2) && (MCAtom[type1].molecule == 2)) && (MCAtom[IMTYPE].numb > 1) )
         {
          //   if ( (MCType[atom0] == IMTYPE) && (MCType[atom1] == IMTYPE) )
         //   {
         //     if ( (MCAtom[type0].molecule == 2) && (MCAtom[type1].molecule == 2) )
         //   {
        //     cout<<"Le if de GG: MCAtom[type0].molecule MCAtom[type1].molecule"<<MCAtom[type0].molecule<<" "<<MCAtom[type1].molecule<<endl;
             double com_1[3];
             double com_2[3];
             double Eulang_1[3];
             double Eulang_2[3];
             double E_2H2O;
             for (int id=0;id<NDIM;id++)
             {
                  com_1[id] = MCCoords[id][t0];
                  com_2[id] = MCCoords[id][t1];
             }
             int tm0=offset0 + it/RotRatio;
             int tm1=offset1 + it/RotRatio;
             Eulang_1[PHI]=MCAngles[PHI][tm0];
             Eulang_1[CTH]=acos(MCAngles[CTH][tm0]);
             Eulang_1[CHI]=MCAngles[CHI][tm0];
             Eulang_2[PHI]=MCAngles[PHI][tm1];
             Eulang_2[CTH]=acos(MCAngles[CTH][tm1]);
             Eulang_2[CHI]=MCAngles[CHI][tm1];
             caleng_(com_1, com_2, &E_2H2O,
                        Eulang_1, Eulang_2);
             spot_pair += E_2H2O;
          //  }  //
          //  }
         }
//------------- [ATOM - ATOM] ------------------------------- 
         else                // only one atom type    
         if ((type0 == type1) && MCAtom[type0].molecule == 0) // no "cross" densities 
         {
//          bin_1Ddensity (r,type1);    // densities 
            spot_pair += SPot1D(r,type1);    // potential energy
         }
      }  // LOOP OVER TIME SLICES
      spot += spot_pair;
   }     // LOOP OVER ATOM PAIRS

// cout<<"in GetPotDensity"<<" _gr1D[0][80]="<<_gr1D[0][80]<<" _gr1D_sum[0][80]="<<_gr1D_sum[0][80]<<endl;
   return (spot/(double)NumbTimes);
}

double GetKinEnergy(void)
{
#ifdef DEBUG_PIMC
   const char *_proc_=__func__; //  GetKinEnergy() 
#ifdef DEBUG_WORM 
   if (Worm.exists)
   nrerror(_proc_," Only for Z-configurations");
#endif
#endif

   int    numb  = 0;       // atom number counter, grand canonical
   double r2avr = 0.0;     // <r^2> 

   for (int atom=0;atom<NumbAtoms;atom++)
   {  
      numb ++;            // grand canonical only
 
      int type    = MCType[atom];
      int offset0 = NumbTimes*atom;
      int offset1;    
      
      int gatom   = MCAtom[type].offset/NumbTimes;   

      double sum = 0.0;
 
      #pragma omp parallel for reduction(+: sum)
      for (int it=0;it<NumbTimes;it++) 
      {
          int t0  = offset0 + it;
        
          offset1 = offset0;
          if ((MCAtom[type].stat ==  BOSE) && ((it+1) == NumbTimes))          
          offset1 = NumbTimes*(gatom + PIndex[atom-gatom]);

	  int t1  = offset1 + (it+1) % NumbTimes; // = offset1

          for (int dim=0;dim<NDIM;dim++)
	  {
             double dr = MCCoords[dim][t0] - MCCoords[dim][t1];

             if (MINIMAGE)
             dr  -= (BoxSize[dim]*rint(dr/BoxSize[dim]));

             sum += (dr*dr);
          }    
       } // END loop over time slices

       r2avr += (sum/(4.0*MCBeta*MCAtom[type].lambda)); 

   }    // END loop over atoms

#ifdef DEBUG_PIMC
   if (numb != NumbAtoms)   // should be removed for grand canonical calculations            
   nrerror(_proc_,"Wrong number of atoms");
#endif

// r2avr /= (double)numb;

   double kin = (double)NumbTimes*Temperature*(0.5*(double)(NDIM*numb) - r2avr);

   return kin;
}

double GetRotEnergy(void)
{
   int type = IMTYPE; 

   int offset = MCAtom[type].offset;
 
   int atom  = 0;                   // only one molecular impurtiy
   offset   += (NumbTimes*atom);
   int gatom = offset/NumbTimes;    // the same offset for rot and trans degrees

   double srot = 0.0;
   ErotSQ=0.0;
   Erot_termSQ=0.0;

   for (int it0=0;it0<NumbRotTimes;it0++)
   {
      int t0 = offset +  it0;
      int t1 = offset + (it0 + 1) % NumbRotTimes;

      double p0 = 0.0;
      for (int id=0;id<NDIM;id++)
      p0 += (MCCosine[id][t0]*MCCosine[id][t1]);

      if(RotDenType == 0)
      {
         double rdens = SRotDens(p0,type);

         if  (fabs(rdens)>RZERO)               // need to find asymptotic for small rot dens
         srot += (SRotDensDeriv(p0,type)/rdens);
         Erot_termSQ += (SRotDensDeriv(p0,type)/rdens)*(SRotDensDeriv(p0,type)/rdens);
         ErotSQ += SRotDensEsqrt(p0,type)/rdens;
      }
      else if(RotDenType == 1)
      {
         double rho,erot;
         rsline_(&X_Rot,&p0,&MCRotTau,&rho,&erot);
//       srot += erot/(double)NumbRotTimes;
         srot += rho;
      }
   }

   if(RotDenType == 1)
   {
      srot = srot/(double)NumbRotTimes;
      srot = srot/MCRotTau + 1.0/MCRotTau;
   }
   
   return (srot);	      
}

double GetRotE3D(void)
{
   int type = IMTYPE;

   int offset = MCAtom[type].offset;

   double ERot3D=0.0;

   ErotSQ=0.0;
   Erot_termSQ=0.0;

   for(int atom  = 0;atom<MCAtom[type].numb;atom++)                   // multi molecular impurtiy
   {
   offset   += (NumbTimes*atom);
   int gatom = offset/NumbTimes;    // the same offset for rot and trans degrees

   double srot = 0.0;
   double sesq = 0.0;
   double se_termsq=0.0;

   int RNskip;
   if(RotDenType == 0)
   {
      RNskip = 1;
   }
   else if(RotDenType == 1)
   {
      RNskip = RNratio;
   }
   
   #pragma omp parallel for reduction(+: srot,sesq)
   for (int it0=0;it0<NumbRotTimes;it0=it0+RNskip)
   {
      int t0 = offset +  it0;
      int t1 = offset + (it0 + RNskip) % NumbRotTimes;
//    Given the two sets of Euler angles at t0 and t1, Toby calculates srot and sesq

      double rho;
      double erot;
      double esq;
      int istop=0;
      double Eulan1[3];
      double Eulan2[3];
      double Eulrel[3];
      double therel,phirel,chirel;

      Eulan1[0]=MCAngles[PHI][t0];
      Eulan1[1]=acos(MCAngles[CTH][t0]);
      Eulan1[2]=MCAngles[CHI][t0];
      Eulan2[0]=MCAngles[PHI][t1];
      Eulan2[1]=acos(MCAngles[CTH][t1]);
      Eulan2[2]=MCAngles[CHI][t1];

      rotden_(Eulan1,Eulan2,Eulrel,&rho,&erot,&esq,rhoprp,erotpr,erotsq,&istop);
      phirel=Eulrel[0];
      therel=Eulrel[1];
      chirel=Eulrel[2];

      int bin_t    = (int)floor(therel/_delta_theta);
      if ((bin_t<MC_BINST) && (bin_t>=0))
      _relthe_sum[bin_t] +=(double)RNskip;

      int bin_p = (int)floor(phirel/_delta_chi);
      if ((bin_p<MC_BINSC) && (bin_p>=0))
      _relphi_sum[bin_p] +=(double)RNskip;

      int bin_c = (int)floor(chirel/_delta_chi);
      if ((bin_c<MC_BINSC) && (bin_c>=0))
      _relchi_sum[bin_c] +=(double)RNskip;

//    Rattle Shake rotational energy
      if(RotDenType == 1 && RNratio == 1)
      rsrot_(Eulan1,Eulan2,&X_Rot,&Y_Rot,&Z_Rot,&MCRotTau,&RotOdEvn,&RotEoff,&rho,&erot);

//    srot += erot;
      if(RotDenType == 1 && RNratio == 1)
      {
         srot += rho;
      }
      else
      srot += erot;

      sesq += esq;
      se_termsq += erot*erot;

   }

   srot = srot / ((double)(NumbRotTimes/RNskip));
   sesq = sesq / ((double)(NumbRotTimes/RNskip)*(double)(NumbRotTimes/RNskip));
   se_termsq = se_termsq/ ((double)(NumbRotTimes/RNskip)*(double)(NumbRotTimes/RNskip));
   ERot3D += srot;
   ErotSQ += sesq;
   Erot_termSQ += se_termsq;

// Rattle Shake rotational energy
   if(RotDenType == 1 && RNratio == 1)
   {
      ERot3D = ERot3D/(4.0*(MCRotTau/WNO2K)*(MCRotTau/WNO2K)); 
      ERot3D += 0.25*(X_Rot+Y_Rot+Z_Rot) + 1.5/(MCRotTau/WNO2K);

      ERot3D = ERot3D/WNO2K;
   }

   }


   return (ERot3D);
}

/* reactive */
void GetRCF(void)
//
//  rotational correlation function
//
{
   int type = IMTYPE; 

   int offset = MCAtom[type].offset; // the same offset for rot and trans coordinates
 
   int atom  = 0;                    // only one molecular impurtiy
   offset   += (NumbTimes*atom);
   int gatom = offset/NumbTimes;

   for (int it0=0;it0<NumbRotTimes;it0++)    
   {
      int t0 = offset +  it0;

      for (int itc=0;itc<NumbRotTimes;itc++)  // offsets
      {	
          int tc = offset + (it0 + itc) % NumbRotTimes;

          double p0 = 0.0;
          for (int id=0;id<NDIM;id++)
          p0 += (MCCosine[id][t0]*MCCosine[id][tc]);

         _rcf     [0][itc] += p0;   // block average
         _rcf_sum [0][itc] += p0;   // total average  

	  for (int in=1;in<NUMB_RCF;in++)       // rcf[0][] should be the same as rcf[1][]
          {
              double pleg = 1.0;
 
	      if (p0<PLONE)
//	      pleg = gsl_sf_legendre_Pl(in,p0); // inefficient 	  

             _rcf     [in][itc] += pleg;     
             _rcf_sum [in][itc] += pleg;     
	  }
      } // END offsets	
   }    // END average over the time origin 
}

void SaveRCF(const char fname [], double acount, int mode)
//  save rotational correlation functions
//  
//  mode:  MC_TOTAL - accumulated averages
//  mode:  MC_BLOCK - block averages
//
{
  fstream fid;
  string  frcf;

  frcf  = fname;

  if (mode == MC_TOTAL)    // accumulated averages
  frcf += IO_SUM; 

  frcf += IO_EXT_RCF;

  fid.open(frcf.c_str(),ios::out); io_setout(fid);

  double norm = acount*(double)NumbRotTimes;

  double ** rcf_save;

  rcf_save = _rcf;
  if (mode == MC_TOTAL)    // accumulated averages
  rcf_save = _rcf_sum;

  for (int it=0;it<=NumbRotTimes;it++)    // save <n(tau)n(0)>
  {	  
     fid << setw(IO_WIDTH) << (double)it*MCRotTau << BLANK; 
     fid << setw(IO_WIDTH) << rcf_save[0][it % NumbRotTimes]/norm << BLANK;
 
     fid << endl;
  }

  fid << endl;  // gnuplot index : at list two blank lines
  fid << endl;
  fid << COMMENTS << endl;

  for (int it=0;it<=NumbRotTimes;it++)              // save <Pl(nn)>
  {	  
     fid << setw(IO_WIDTH) << (double)it*MCRotTau << BLANK; 
 
     for (int ip=1;ip<NUMB_RCF;ip++) 
     fid << setw(IO_WIDTH) << rcf_save[ip][it % NumbRotTimes]/norm << BLANK; 
    
     fid << endl; 
  }

  fid.close();
}
/* reactive */

double GetConfPoten_Densities(void)
// should be compatible with ConfPot() from mc_piqmc.cc
{
   const char *_proc_=__func__; //  GetPotEnergy_Densities()  

   if (Worm.exists)
   nrerror(_proc_," Only for Z-configurations");

   double spot = 0.0;  		
 
   for (int atom=0;atom<NumbAtoms;atom++) 	   
   {
      int offset = NumbTimes*atom;
      int type   = MCType[atom];

      for (int it=0;it<NumbTimes;it++) 	    
      { 
//       bool wline = true;                  // skip if the time slice between ira and masha

//       if (WORM && Worm.exists && (Worm.type == type))  
//       wline = WorldLine((atom1-MCAtom[type1].offset/NumbTimes), it);
          
//       if (wline)
         {
         double r2 = 0.0;  		
         for (int id=0;id<NDIM;id++)
         r2   += (MCCoords[id][offset+it]*MCCoords[id][offset+it]);

         spot += (MCAtom[type].mass*r2);

         bin_1Ddensity (sqrt(r2),type);    // densities 

         }
      }
   }

   return (0.5*HOMEGA*HOMEGA*spot/(double)NumbTimes);
}

void bin_2Ddensity(double r, double cost, int dtype)
// (r, cost)  ==  (radius, cos(theta))
//  dtype     ==   density type
{
   int bin_r = (int)floor((r-_min_radius)/_delta_radius);

   if ((bin_r<MC_BINSR) && (bin_r>=0))
   { 
      double theta = acos(cost);
      int bin_t    = (int)floor(theta/_delta_theta);

      if ((bin_t<MC_BINST) && (bin_t>=0))
      {
         _gr2D[dtype][bin_r][bin_t] += 1.0;  // block average
         _gr2D_sum[dtype][bin_r][bin_t] += 1.0; // total average //added by Hui Li
      }
   }
}

void bin_3Ddensity(double r, double theta, double chi, int dtype)
// dtype == density type
{

   int bin_r = (int)floor((r-_min_radius)/_delta_radius);
   if ((bin_r<MC_BINSR) && (bin_r>=0))
   {
      int bin_t    = (int)floor(theta/_delta_theta);
      if ((bin_t<MC_BINST) && (bin_t>=0))
      {
         int bin_c = (int)floor(chi/_delta_chi);
         if((bin_c < MC_BINSC) && (bin_c >= 0))
         {
           int ijk = ((bin_r * MC_BINST) + bin_t) * MC_BINSC + bin_c;
           _gr3D[dtype][ijk] += 1.0; // block average
           _gr3D_sum[dtype][ijk] += 1.0; // total average
         }
      }

   }

}

void bin_1Ddensity(double r,int dtype)
//  r         ==   radius
//  dtype     ==   density type
{
   int bin_r = (int)floor((r-_min_radius)/_delta_radius);

   if ((bin_r<MC_BINSR) && (bin_r>=0))
   {
     _gr1D[dtype][bin_r] += 1.0;
     _gr1D_sum[dtype][bin_r] += 1.0;
   }
}

void SaveGraSum(const char fname [], double acount)
// accumulate sum for inter-atomic distribution.  should be similar to the pair distribution in SaveDensities1D
{
  fstream fid;
  string fdens;

  fdens  = fname;
  fdens += IO_SUM;
  fdens += IO_EXT_GRA;

  fid.open(fdens.c_str(),ios::out); io_setout(fid);

  double norma = _delta_radius*acount*(double)(NumbTimes);

  double r;

  for (int ir=0;ir<MC_BINSR;ir++)
  {
     r   =  (double)ir*_delta_radius;
     r  +=  (0.5*_delta_radius);
     r  +=  _min_radius;

     fid<<setw(IO_WIDTH)<<(r*Units.length)<<BLANK;

     for (int id=0;id<NUMB_DENS1D;id++)
     {
//      double nfact = norma*(double)(MCAtom[id].numb*(MCAtom[id].numb-1));
//      the following scaling is to let the gra_sum to be normalized to one by integrating over dr, without any jacobian factor
//      norma = norma * (MCAtom[id].numb*(MCAtom[id].numb-1))/2.0;
        fid <<setw(IO_WIDTH)<<_gr1D_sum[id][ir]/(norma*(MCAtom[id].numb*(MCAtom[id].numb-1))/2.0)<<BLANK;   // gra_sum
     }

     fid<<endl;
  }
//cout<<"norma="<<norma<<" _gr1D_sum[0][80]="<<_gr1D_sum[0][80]<<" numb="<<<<endl;

  fid.close();
}

void SaveDensities1D(const char fname [], double acount)
// the density type corresponds to the atom type
{
  fstream fid;
  string fdens;

// ------- 1D radial (pair) distribution functions -----------------------------

  fdens  = fname;
  fdens += IO_EXT_GRA;

  fid.open(fdens.c_str(),ios::out); io_setout(fid);

  double volume = 1.0;
  for (int id=0;id<NDIM;id++)
  volume *= BoxSize[id];

  double norm0  = 2.0*M_PI*_delta_radius*acount  // normalization factor for 
                *(double)(NumbTimes)/volume;     // radial distribution functions 

  double r,r2;

  for (int ir=0;ir<MC_BINSR;ir++) // normalization
  {	  
     r   =  (double)ir*_delta_radius;
     r  +=  (0.5*_delta_radius);
     r  +=  _min_radius;
     r2  =  r*r;
   
     fid<<setw(IO_WIDTH)<<(r*Units.length)<<BLANK; 
 
     for (int id=0;id<NUMB_DENS1D;id++)
     { 
        double nfact = norm0*(double)(MCAtom[id].numb*(MCAtom[id].numb-1));
        fid <<setw(IO_WIDTH)<<_gr1D[id][ir]/(r2*nfact)<<BLANK;   // gra
     } 

     fid<<endl;
  }

  fid.close();

// CONVERT 2D TO 1D
  
  if (IMPURITY && MCAtom[IMTYPE].molecule == 1)  // ? separate from 1D case to avoid long if{}  ?
  {
// ------  1D radial distributions around an impurity --------------------------

  fdens  = fname;
  fdens += IO_EXT_GRI;

  fid.open(fdens.c_str(),ios::out); io_setout(fid);

  double  norm1 = _delta_radius*(double)NumbTimes     // the normalization for the RADIAL
                  *acount;                            // distribution around impurity

  norm1 *= (4.0*M_PI*pow(Units.length,(double)NDIM)); // 3D only

  for (int ir=0;ir<MC_BINSR;ir++) // radial distribution
  {
     r   =  (double)ir*_delta_radius;
     r  +=  (0.5*_delta_radius);
     r  +=  _min_radius;
     r2  =   r*r;
   
     fid << setw(IO_WIDTH) << (r*Units.length) << BLANK; 

     for (int id=0;id<NUMB_DENS1D;id++) // rho(r,cost) -> rho(r) and rho(t) convert 
     {
        double densr = 0.0;          
        for (int it=0;it<MC_BINST;it++)
        densr += _gr2D[id][ir][it]; 
  
        fid << setw(IO_WIDTH) << densr/(norm1*r2) << BLANK;  
     }
      
     fid << endl;
   } 

   fid.close();

//fid << endl;  // gnuplot index : at list two blank lines
//fid << endl;

// ------  1D angular distributions around an impurity --------------------------

   fdens  = fname;
   fdens += IO_EXT_GRT;

   fid.open(fdens.c_str(),ios::out); io_setout(fid);

   double norm2 = (double)NumbTimes*acount*_delta_theta;

   double theta;

   for (int it=0;it<MC_BINST;it++)  // angular distribution
   {
      theta  = (double)it*_delta_theta;
      theta += (0.5*_delta_theta);

      fid << setw(IO_WIDTH) << (theta*180.0/M_PI) << BLANK;

      for (int id=0;id<NUMB_DENS1D;id++) // rho(r,cost) -> rho(r) and rho(t) convert 
      {
         double denst = 0.0;
         for (int ir=0;ir<MC_BINSR;ir++)  // total contribution
         denst += _gr2D[id][ir][it];

         fid << setw(IO_WIDTH) << denst/(norm2*(double)MCAtom[id].numb) << BLANK; 
      }

      fid << endl;
   }
 
   fid.close();
   } // END if (IMPURITY) convert 2D into 1D
}

void SaveRho1D(const char fname [], double acount, int mode)
// the density type corresponds to the atom type
// this subroutine is designed for nonlinear rotor only and should not be called for linear rotor under
// ALL circumstances
{
  fstream fid;
  string fdens;

  double ** dens_save;

  dens_save = _gr3D;
  if(mode == MC_TOTAL)
  dens_save = _gr3D_sum;


// convert 3D density to 1D and save
  if (IMPURITY && MCAtom[IMTYPE].molecule == 2)  //
  {
// ------  1D radial distributions around an impurity --------------------------

  fdens  = fname;

  if(mode == MC_TOTAL) // accumulated averages
  fdens += IO_SUM;

  fdens += IO_EXT_GRI;

  fid.open(fdens.c_str(),ios::out); io_setout(fid);

  double  norm1 = _delta_radius*(double)NumbTimes     // the normalization for the RADIAL
                  *acount;                            // distribution around impurity

  double r,r2;

  for (int ir=0;ir<MC_BINSR;ir++) // radial distribution
  {
     r   =  (double)ir*_delta_radius;
     r  +=  (0.5*_delta_radius);
     r  +=  _min_radius;
     r2  =   r*r;

     fid << setw(IO_WIDTH) << (r*Units.length) << BLANK;

     for (int id=0;id<NUMB_DENS3D;id++) // rho(r,theta,chi) -> rho(r), rho(theta) and rho(chi) convert 
     {
        double densr = 0.0;
        for (int it=0;it<MC_BINST;it++)
        for (int ic=0;ic<MC_BINSC;ic++)
        {
           int ijk = (ir*MC_BINST + it)*MC_BINSC + ic;
           densr += dens_save[id][ijk];
        }

//      fid << setw(IO_WIDTH) << densr/(norm1*r2) << BLANK;  
        fid << setw(IO_WIDTH) << densr/(norm1) << BLANK;

     }

     fid << endl;
   }

   fid.close();

//fid << endl;  // gnuplot index : at list two blank lines
//fid << endl;

// ------  1D theta distributions around an impurity --------------------------

   fdens  = fname;

   if(mode == MC_TOTAL) // accumulated averages
   fdens += IO_SUM;

   fdens += IO_EXT_GRT;

   fid.open(fdens.c_str(),ios::out); io_setout(fid);

   double norm2 = (double)NumbTimes*acount*_delta_theta*(180.0/M_PI); // the last factor converts density to per degree

   double theta;

   for (int it=0;it<MC_BINST;it++)  // angular distribution
   {
      theta  = (double)it*_delta_theta;
      theta += (0.5*_delta_theta);

      fid << setw(IO_WIDTH) << (theta*180.0/M_PI) << BLANK;

      for (int id=0;id<NUMB_DENS3D;id++) // rho(r,cost) -> rho(r) and rho(t) convert 
      {
         double denst = 0.0;
         for (int ir=0;ir<MC_BINSR;ir++)  // total contribution
         for (int ic=0;ic<MC_BINSC;ic++)
         {
           int ijk = (ir*MC_BINST + it)*MC_BINSC + ic;
           denst += dens_save[id][ijk];
         }

         fid << setw(IO_WIDTH) << denst/(norm2*(double)MCAtom[id].numb) << BLANK;
      }

      fid << endl;
   }

   fid.close();

// ------  1D chi distributions around an impurity --------------------------
   fdens  = fname;

   if(mode == MC_TOTAL) // accumulated averages
   fdens += IO_SUM;

   fdens += IO_EXT_GRC;

   fid.open(fdens.c_str(),ios::out); io_setout(fid);

   double norm4 = (double)NumbTimes*acount*_delta_chi*(180./M_PI); // the last factor is to convert the density in the chi element of degree.

   double chi;

   for (int ic=0;ic<MC_BINSC;ic++)
   {
      chi  = (double)ic*_delta_chi;
      chi += (0.5*_delta_chi);

      fid << setw(IO_WIDTH) << (chi*180.0/M_PI) << BLANK;

      for (int id=0;id<NUMB_DENS3D;id++)
      {
          double densc = 0.0;
          for (int ir=0;ir<MC_BINSR;ir++)
          for (int it=0;it<MC_BINST;it++)
          {
           int ijk = (ir*MC_BINST + it)*MC_BINSC + ic;
           densc += dens_save[id][ijk];
          }

         fid << setw(IO_WIDTH) << densc/(norm4*(double)MCAtom[id].numb) << BLANK;
      }

      fid << endl;
   }

   fid.close();

// ------  Punch out relative euler angles between the adjacent rotational imaginary time slices  --------------------------
      if(mode == MC_TOTAL)
      {
         fdens = fname;
         fdens += IO_SUM;
         fdens += IO_EXT_REP;

         fid.open(fdens.c_str(),ios::out); io_setout(fid);

         double norm5 = (double)NumbRotTimes*acount*_delta_chi*(180./M_PI);

         double phirel;

         for (int ip=0;ip<MC_BINSC;ip++)
         {
            phirel = (double)ip*_delta_chi;
            phirel += (0.5*_delta_chi);

            fid << setw(IO_WIDTH) << (phirel*180.0/M_PI) << BLANK;

            fid << setw(IO_WIDTH) << _relphi_sum[ip]/norm5 << endl;

         }

         fid.close();

         fdens = fname;
         fdens += IO_SUM;
         fdens += IO_EXT_REC;

         fid.open(fdens.c_str(),ios::out); io_setout(fid);

         double chirel;

         for (int ic=0;ic<MC_BINSC;ic++)
         {
            chirel = (double)ic*_delta_chi;
            chirel += (0.5*_delta_chi);

            fid << setw(IO_WIDTH) << (chirel*180.0/M_PI) << BLANK;

            fid << setw(IO_WIDTH) << _relchi_sum[ic]/norm5 << endl;

         }

         fid.close();

         fdens = fname;
         fdens += IO_SUM;
         fdens += IO_EXT_RET;

         fid.open(fdens.c_str(),ios::out); io_setout(fid);

         double norm6 = (double)NumbRotTimes*acount*_delta_theta*(180.0/M_PI);

         double therel;

         for (int it=0;it<MC_BINST;it++)
         {
            therel = (double)it*_delta_theta;
            therel += (0.5*_delta_theta);

            fid << setw(IO_WIDTH) << (therel*180.0/M_PI) << BLANK;

            fid << setw(IO_WIDTH) << _relthe_sum[it]/norm6 << endl;

         }

         fid.close();

      }

   } // END if (IMPURITY) convert 3D into 1D
}

// added by Hui Li
void SaveDensities2D(const char fname [], double acount, int mode)
// the density type corresponds to the atom type
//
//  mode:  MC_TOTAL - accumulated averages
//  mode:  MC_BLOCK - block averages
//
{
  fstream fid;
  string fdens;
         
  if (IMPURITY)  // ? separate from 1D case to avoid long if{}  ?
  {
// ------  1D radial distributions around an impurity --------------------------

  fdens  = fname;

  if (mode == MC_TOTAL)    // accumulated averages
  fdens += IO_SUM;

  fdens += IO_EXT_DENS2D;

  fid.open(fdens.c_str(),ios::out); io_setout(fid);

  double  norm3 = _delta_radius*_delta_theta*(double)NumbTimes     // the normalization for density 
                  *acount;                                        // distribution around impurity

  double *** dens_save;
  
  dens_save = _gr2D;
  if (mode == MC_TOTAL)    // accumulated averages density
  dens_save = _gr2D_sum;


  double theta,r,r2;

  for (int it=0;it<MC_BINST;it++)  // angular distribution
   {
      theta  = (double)it*_delta_theta;
      theta += (0.5*_delta_theta);

      for (int ir=0;ir<MC_BINSR;ir++) // radial distribution
       {
         r   =  (double)ir*_delta_radius;
         r  +=  (0.5*_delta_radius);
         r  +=  _min_radius;
         r2  =   r*r;
    
         for (int id=0;id<NUMB_DENS2D;id++) // 
          {
           fid << setw(IO_WIDTH) << (theta*180.0/M_PI) << BLANK;
           fid << setw(IO_WIDTH) << (r*Units.length) << BLANK;
           fid << setw(IO_WIDTH) << dens_save[id][ir][it]/(norm3) << BLANK;
          }

         fid << endl;
        }
   }
   fid.close();
  } //endif 
 } //end subroutine

// added by Toby Zeng
void SaveDensities3D(const char fname [], double acount, int mode)
// the density type corresponds to the atom type
//
//  mode:  MC_TOTAL - accumulated averages
//  mode:  MC_BLOCK - block averages
//
{
  fstream fid;
  string fdens;

  if (IMPURITY)  // ? separate from 1D case to avoid long if{}  ?
  {
// ------  1D radial distributions around an impurity --------------------------

  fdens  = fname;

  if (mode == MC_TOTAL)    // accumulated averages
  fdens += IO_SUM;

  fdens += IO_EXT_DENS3D;

  fid.open(fdens.c_str(),ios::out); io_setout(fid);

  double  norm5 = _delta_radius*_delta_theta*_delta_chi*(double)NumbTimes     // the normalization for density 
                  *acount;                                        // distribution around impurity

  double ** dens_save;

  dens_save = _gr3D;
  if (mode == MC_TOTAL)    // accumulated averages density
  dens_save = _gr3D_sum;


  double theta,r,r2,chi;

  for (int ir=0;ir<MC_BINSR;ir++)
  {
      r   =  (double)ir*_delta_radius;
      r  +=  (0.5*_delta_radius);
      r  +=  _min_radius;

      for (int it=0;it<MC_BINST;it++)
      {
         theta  = (double)it*_delta_theta;
         theta += (0.5*_delta_theta);

         for (int ic=0;ic<MC_BINSC;ic++)
         {
            chi  = (double)ic*_delta_chi;
            chi += (0.5*_delta_chi);

            for (int id=0;id<NUMB_DENS3D;id++)
            {
               fid << setw(IO_WIDTH) << (r*Units.length) << BLANK;
               fid << setw(IO_WIDTH) << (theta*180.0/M_PI) << BLANK;
               fid << setw(IO_WIDTH) << (chi*180.0/M_PI) << BLANK;
               int ijk = (ir*MC_BINST + it)*MC_BINSC + ic;
               fid << setw(IO_WIDTH) << dens_save[id][ijk]/(norm5*r*r*sin(theta)) << BLANK;
            }
            fid<<endl;
         }

      }
  }

/*
  for (int it=0;it<MC_BINST;it++)  // angular distribution
   {
      theta  = (double)it*_delta_theta;
      theta += (0.5*_delta_theta);

      for (int ir=0;ir<MC_BINSR;ir++) // radial distribution
       {
         r   =  (double)ir*_delta_radius;
         r  +=  (0.5*_delta_radius);
         r  +=  _min_radius;
         r2  =   r*r;

         for (int id=0;id<NUMB_DENS2D;id++) // 
          {
           fid << setw(IO_WIDTH) << (theta*180.0/M_PI) << BLANK;
           fid << setw(IO_WIDTH) << (r*Units.length) << BLANK;
           fid << setw(IO_WIDTH) << dens_save[id][ir][it]/(norm5) << BLANK;
          }

         fid << endl;
        }
   }
*/
   fid.close();
  } //endif 
 } //end subroutine

// added by Toby Zeng
void IOxyzAng(int tstatus, const char file_name[])
{
   const char *_proc_=__func__;    // "IOxyz"; 

   string fdens;

//---------------- Open  ------------

   ios::openmode mode;

   switch (tstatus)
   {
      case IOWrite: mode = ios::out;  break;
      case IORead : mode = ios::in;   break;
      default     :
      nrerror (_proc_,IO_ERR_WMODE);  break;
   }

   fdens = file_name;

   if(IOWrite)
   fdens += IO_EXT_XYZ;

   fstream fid(fdens.c_str(),mode);

   if (!fid.good())
   _io_error(_proc_,IO_ERR_FOPEN,file_name);

   io_setout(fid);

//---------------- Read/Write ------------

// stringstream stype;
   string       sbuff,stype;

   int offset;

   int type = 0;
   int atom = 0;  // first atom # will be 1, NOT 0 
   switch (tstatus)
   {
      case IOWrite:
         fid<<MaxnTimes<<" ";                    // total number of "atoms" 
//       permutation table
         if(BOSONS)
         {
            for (int atom=0; atom<MCAtom[BSTYPE].numb;atom++)
            {
             fid <<"  "<<PIndex[atom]<<BLANK;
            }
         }
         fid << endl;

         fid<<COMMENTS<<BLANK<<IO_COM_XYZ<<endl;  // comments

         for (int it=0;it<MaxnTimes;it++)
         {
            if (it==MCAtom[type+1].offset) {type++; atom=0;}    // new atom type
            if ((it-MCAtom[type].offset)%NumbTimes==0) atom++;  // new atom

            fid<<MCAtom[type].type<<atom;      // atom label
//          fid<<setw(5)<<stype<<BLANK;

            for (int id=0;id<NDIM;id++)
            {
            fid<<setw(IO_WIDTH)<<MCCoords[id][it]<<BLANK;
            fid<<setw(IO_WIDTH)<<MCAngles[id][it]<<BLANK;
//          Toby replaces the above line by
//          fid<<setw(IO_WIDTH)<<MCAngles[id][it]<<BLANK;
//          by doing that, Toby stores the three Euler angles, not the unit vector of the two angles orientation
            }
            fid<<endl;
         }

         break;

      case IORead:
         fid>>MaxnTimes;
         getline(fid,sbuff);    // skip a comment line

         for (int type=0;type<NumbTypes;type++)
         for (int atom=0;atom<MCAtom[type].numb;atom++)
         for (int it=0;it<NumbTimes;it++)
         {
            fid>>sbuff;         // skip an atom type

            offset=MCAtom[type].offset+NumbTimes*atom;
            for (int id=0;id<NDIM;id++)
            {
             fid>>MCCoords[id][offset+it];
             fid>>MCCosine[id][offset+it];
//           Toby replaces the above line by
//           fid>>MCAngles[id][offset+it];
//           by doing that, Toby reads the three Euler angles, not the unit vector of the two angle orientation
            }
         }

         break;

      default :
         nrerror (_proc_,IO_ERR_WMODE);
         break;
   }

   fid.close();
}

// added by Toby Zeng
void SaveRhoThetaChi(const char fname [], double acount, int mode)
// the density type corresponds to the atom type
//
//  mode:  MC_TOTAL - accumulated averages
//  mode:  MC_BLOCK - block averages
//
{
  fstream fid;
  string fdens;

  fdens  = fname;

  if (mode == MC_TOTAL)    // accumulated averages
  fdens += IO_SUM;

  fdens += IO_EXT_GTC;

  fid.open(fdens.c_str(),ios::out); io_setout(fid);

  double  norm6 = _delta_theta*_delta_chi*(double)NumbTimes     // the normalization for density 
                  *acount*(180.0*180.0/(M_PI*M_PI));            // distribution around impurity

  double ** dens_save;

  dens_save = _gr3D;
  if (mode == MC_TOTAL)    // accumulated averages density
  dens_save = _gr3D_sum;


  double theta,r,r2,chi;

  for (int it=0;it<MC_BINST;it++)
  {
     theta  = (double)it*_delta_theta;
     theta += (0.5*_delta_theta);

     for (int ic=0;ic<MC_BINSC;ic++)
     {
        chi  = (double)ic*_delta_chi;
        chi += (0.5*_delta_chi);

        fid << setw(IO_WIDTH) << (theta*180.0/M_PI) << BLANK;
        fid << setw(IO_WIDTH) << (chi*180.0/M_PI) << BLANK;

        for (int id=0;id<NUMB_DENS3D;id++)
        {

           double denstc = 0.0;

           for (int ir=0;ir<MC_BINSR;ir++)
           {
              int ijk = (ir*MC_BINST + it)*MC_BINSC + ic;
              denstc = denstc + dens_save[id][ijk];
//            fid << setw(IO_WIDTH) << dens_save[id][ijk]/(norm6) << BLANK;
           }
           fid << setw(IO_WIDTH) << denstc/(norm6*(double)MCAtom[id].numb) << BLANK;
        }
        fid<<endl;
     }

     fid<<endl;

  }

   fid.close();
 } //end subroutine

void GetExchangeLength(void)
{
   for (int atom=0;atom<MCAtom[BSTYPE].numb;atom++)
  _pflags[atom] = 0;

   for (int atom=0;atom<MCAtom[BSTYPE].numb;atom++)
   if (_pflags[atom] == 0)
   { 
//    _pflags[patom] = 1;           // do not need if try all atoms in order
 
       int clenght = 0;             // start a new cycle
       int patom   = PIndex[atom]; 
	
       while (patom != atom)
       {
         _pflags[patom] = 1; 
          patom = PIndex[patom];
          clenght++;
       }

      _ploops[clenght] += 1.0;	 
   }
}

void SaveExchangeLength (const char fname [], double acount, long int blocknumb)
{
  const char *_proc_=__func__;    //  SaveExchangeLength() 
 
//------  open file ---------------

  fstream fid;
  string  fperm;

  fperm  = fname;
  fperm += IO_EXT_PRL;

  fid.open(fperm.c_str(),ios::app | ios::out); io_setout(fid);

  if (!fid.is_open())
 _io_error(_proc_,IO_ERR_FOPEN,fperm.c_str());

//--------------------------------

  fid << setw(IO_WIDTH_BLOCK) << blocknumb << BLANK;     // block number

  double excited = 0.0;  
  double ground  = 0.0;   

  for (int clength=0;clength<MCAtom[BSTYPE].numb;clength++)
  {
//   check <norm> below
     double norm  = (double)(clength+1)/(acount*(double)MCAtom[BSTYPE].numb);

     if (clength<=GSLOOP_MAX)  
     excited  += (_ploops[clength]*norm);
     else
     ground   += (_ploops[clength]*norm);
  }

  fid << setw(IO_WIDTH) << ground  << BLANK;   
  fid << setw(IO_WIDTH) << excited << BLANK;  
 
  fid << setw(IO_WIDTH) << (ground+excited) << BLANK; // norm check

// ----------------------------------------------

  for (int clength=0;clength<MCAtom[BSTYPE].numb;clength++)
  {
//   check <norm> above 
     double norm  = (double)(clength+1)/(acount*(double)MCAtom[BSTYPE].numb);

     fid << setw(IO_WIDTH) << _ploops[clength]*norm << BLANK;   
  }
  fid << endl;

//added by Hui Li, test whether it is changed
  for (int atom=0; atom<MCAtom[BSTYPE].numb;atom++)
  {
   fid << setw(IO_WIDTH)<<PIndex[atom]<<BLANK;  
  }

// print PrintXYZprl in the prl file, to specify whether the instantaneous XYZ and PRL are recorded
   fid << PrintXYZprl;

   fid << endl;
// end of test

  fid.close();
}

void GetAreaEstimators(void)
{
   double dr0[NDIM];
   double dr1[NDIM];

   double n_perp[NDIM];
   double n_parl[NDIM];

   double area[NDIM];

   double rn0[NDIM];     // r x n cross product
   double rn1[NDIM];

   double area_perp;
   double area_parl;

   double inert_perp;
   double inert_parl;

// define the ref point for the area estimator
// (i)   instanteneous COM
// (ii)  instanteneous COM of a dopant molecule

// uncomment one section below

/*
//(i)

// double tmass  = 0.0;   // define total mass globally (canonical only)?

// define mass globally change the order of loops

// for (int type=0;type<NumbTypes;type++)
// tmass += (MCAtom[type].mass * (double)MCAtom[type].numb);

   for (int dim=0;dim<NDIM;dim++) // instantaneous center of mass
   for (int it=0;it<NumbTimes;it++)
   {	
      double      tmass  = 0.0;   // define total mass globally (canonical only)?
      newcoords[dim][it] = 0.0;   // center of mass
 
      for (int atom=0;atom<NumbAtoms;atom++)
      {
         int    type = MCType[atom];
         double mass = MCAtom[type].mass;
 
         newcoords[dim][it] += (mass*MCCoords[dim][atom*NumbTimes + it]);	 
         tmass              +=  mass;
      }

      newcoords[dim][it] /= tmass; 
   }
*/

//(ii) only one molecule in the system
   for (int dim=0;dim<NDIM;dim++) // instantaneous center of mass
   for (int it=0;it<NumbTimes;it++)
   newcoords[dim][it] = MCCoords[dim][MCAtom[IMTYPE].offset + it];	 
   int offset = MCAtom[BSTYPE].offset;

   area_perp  = 0.0;
   area_parl  = 0.0;

   inert_perp = 0.0;
   inert_parl = 0.0;

   for (int atom=0;atom<MCAtom[BSTYPE].numb;atom++)
   for (int it0=0;it0<NumbTimes;it0++)
   {
      int it1 = (it0 + 1) % NumbTimes;

      int pt0 = offset + NumbTimes*atom;     // offset only
      int pt1 = pt0;
 
      if (it1!= (it0 + 1))          
      pt1 = offset + NumbTimes*PIndex[atom]; // offset only  
 
      pt0 += it0;
      pt1 += it1;

      for (int dim=0;dim<NDIM;dim++)         // COM adjustment 
      {
         dr0[dim] = MCCoords[dim][pt0] - newcoords[dim][it0];
         dr1[dim] = MCCoords[dim][pt1] - newcoords[dim][it1];
      }
// ------------- orientations ----------------------------

      int it_rot = it0/RotRatio;

      n_parl[AXIS_X] = MCCosine[AXIS_X][MCAtom[IMTYPE].offset + it_rot];
      n_parl[AXIS_Y] = MCCosine[AXIS_Y][MCAtom[IMTYPE].offset + it_rot];
      n_parl[AXIS_Z] = MCCosine[AXIS_Z][MCAtom[IMTYPE].offset + it_rot];
/*
      n_parl[AXIS_X] = 0.0;
      n_parl[AXIS_Y] = 0.0;
      n_parl[AXIS_Z] = 1.0;
*/
      double static zero = 10e-4; 

      double tg = 0.0;
      double st = 1.0;

      if (fabs(n_parl[AXIS_X]) > zero) // treat separetely n_x=n_y=0 limit ?
      {
         tg = n_parl[AXIS_Y]/n_parl[AXIS_X]; 
         st = sqrt(1.0 + tg*tg);
      }
 
      n_perp[AXIS_X] =  tg /st;
      n_perp[AXIS_Y] = -1.0/st;
      n_perp[AXIS_Z] =  0.0;

/* 
      n_perp[AXIS_X] =  0.0;
      n_perp[AXIS_Y] = -1.0;
      n_perp[AXIS_Z] =  0.0;
*/

//-------------------------------------------------

      area[AXIS_X]  = 0.5*(dr0[AXIS_Y]*dr1[AXIS_Z] - dr0[AXIS_Z]*dr1[AXIS_Y]);
      area[AXIS_Y]  = 0.5*(dr0[AXIS_Z]*dr1[AXIS_X] - dr0[AXIS_X]*dr1[AXIS_Z]);
      area[AXIS_Z]  = 0.5*(dr0[AXIS_X]*dr1[AXIS_Y] - dr0[AXIS_Y]*dr1[AXIS_X]);

      for (int dim=0;dim<NDIM;dim++) // projections
      {
         area_perp += (n_perp[dim]*area[dim]);  
         area_parl += (n_parl[dim]*area[dim]);   	 
      } 

// -- moment of inertia  -----------------------

      rn0[AXIS_X] = n_perp[AXIS_Y]*dr0[AXIS_Z] - n_perp[AXIS_Z]*dr0[AXIS_Y];
      rn0[AXIS_Y] = n_perp[AXIS_Z]*dr0[AXIS_X] - n_perp[AXIS_X]*dr0[AXIS_Z];
      rn0[AXIS_Z] = n_perp[AXIS_X]*dr0[AXIS_Y] - n_perp[AXIS_Y]*dr0[AXIS_X];

      rn1[AXIS_X] = n_perp[AXIS_Y]*dr1[AXIS_Z] - n_perp[AXIS_Z]*dr1[AXIS_Y];
      rn1[AXIS_Y] = n_perp[AXIS_Z]*dr1[AXIS_X] - n_perp[AXIS_X]*dr1[AXIS_Z];
      rn1[AXIS_Z] = n_perp[AXIS_X]*dr1[AXIS_Y] - n_perp[AXIS_Y]*dr1[AXIS_X];

      for (int dim=0;dim<NDIM;dim++) // the perpendicular component
      inert_perp += (rn0[dim]*rn1[dim]);
 
      rn0[AXIS_X] = n_parl[AXIS_Y]*dr0[AXIS_Z] - n_parl[AXIS_Z]*dr0[AXIS_Y];
      rn0[AXIS_Y] = n_parl[AXIS_Z]*dr0[AXIS_X] - n_parl[AXIS_X]*dr0[AXIS_Z];
      rn0[AXIS_Z] = n_parl[AXIS_X]*dr0[AXIS_Y] - n_parl[AXIS_Y]*dr0[AXIS_X];

      rn1[AXIS_X] = n_parl[AXIS_Y]*dr1[AXIS_Z] - n_parl[AXIS_Z]*dr1[AXIS_Y];
      rn1[AXIS_Y] = n_parl[AXIS_Z]*dr1[AXIS_X] - n_parl[AXIS_X]*dr1[AXIS_Z];
      rn1[AXIS_Z] = n_parl[AXIS_X]*dr1[AXIS_Y] - n_parl[AXIS_Y]*dr1[AXIS_X];

      for (int dim=0;dim<NDIM;dim++)  // the parallel component
      inert_parl += (rn0[dim]*rn1[dim]);
   }

  _areas[PERP] += area_perp;   
  _areas[PARL] += area_parl; 
 
  _area2[PERP] += (area_perp*area_perp);   
  _area2[PARL] += (area_parl*area_parl);  

  _inert[PERP] += (inert_perp/(double)NumbTimes);
  _inert[PARL] += (inert_parl/(double)NumbTimes);
}

void GetAreaEstim3D(int iframe)
{

// iframe specifies the frame on which the superfluid response is projected 0: Space-fixed frame; 1: dopant-fixed frame

   double dr0[NDIM];     // the position vector of a bead with respect to the COM of dopant
   double dr1[NDIM];     // the position vector of the bead one slice after

   double inert3D[NDIM*NDIM];

   double rn0[NDIM];     // dr0 x n cross product
   double rn1[NDIM];     // dr1 x n cross product

   double area[NDIM];    // the vector area of dr0 cross dr1

   double area_proj[NDIM]; // projection of the area vector onto the three principle axes of the non-linear dopant

// define the ref point for the area estimator
// (i)   instanteneous COM
// (ii)  instanteneous COM of a dopant molecule

// uncomment one section below

//(i)

// double tmass  = 0.0;   // define total mass globally (canonical only)?

// define mass globally change the order of loops

// for (int type=0;type<NumbTypes;type++)
// tmass += (MCAtom[type].mass * (double)MCAtom[type].numb);
   if(iframe == 0)
   {
///*
   for (int dim=0;dim<NDIM;dim++) // instantaneous center of mass
   for (int it=0;it<NumbTimes;it++)
   {    
      double      tmass  = 0.0;   // define total mass globally (canonical only)?
      newcoords[dim][it] = 0.0;   // center of mass
 
      for (int atom=0;atom<NumbAtoms;atom++)
      {
         int    type = MCType[atom];
         double mass = MCAtom[type].mass;
 
         newcoords[dim][it] += (mass*MCCoords[dim][atom*NumbTimes + it]);        
         tmass              +=  mass;
      }

      newcoords[dim][it] /= tmass; 
   }
//*/
   }
//(ii) only one molecule in the system

   if(iframe == 1)
   {
      for (int dim=0;dim<NDIM;dim++) // instantaneous center of mass
      for (int it=0;it<NumbTimes;it++)
      newcoords[dim][it] = MCCoords[dim][MCAtom[IMTYPE].offset + it];
   }


   int offset = MCAtom[BSTYPE].offset;

   double bmass = MCAtom[BSTYPE].mass;

// clean area_proj and inert3D
   for (int id=0;id<NDIM;id++)
   {
      area_proj[id]=0.0;
      for (int jd=0;jd<NDIM;jd++)
      {
         int ij = id *NDIM + jd;
         inert3D[ij]=0.0;
      }
   }

// private variables that accumulate results in each cpu
   double area_projx=0.0,area_projy=0.0,area_projz=0.0,icxx=0.0,icxy=0.0,icxz=0.0,icyx=0.0,icyy=0.0,icyz=0.0,iczx=0.0,iczy=0.0,iczz=0.0;

// compared to the original GetAreaEstimators, here Toby switches the order of loop over it0 and loop over atom.
// This is to minimize the call of vcord
   #pragma omp parallel for reduction(+: area_projx,area_projy,area_projz,icxx,icxy,icxz,icyx,icyy,icyz,iczx,iczy,iczz) private(area,rn0,rn1,dr0,dr1)
   for (int it0=0;it0<NumbTimes;it0++)
   {

      int it_rot = it0/RotRatio + MCAtom[IMTYPE].offset; // only one molecular impurity

      double RCOM[3];
      double Rpt[3];
      double Eulang[3];
      double vpot3d;
      double radret;
      double theret;
      double chiret;
      double hatx[3];
      double haty[3];
      double hatz[3];
      int    ivcord = 1;

      for (int id=0;id<NDIM;id++)
      {
//       Rpt[id]  = MCCoords[id][it_rot];
//       RCOM[id] = MCCoords[id][it_rot];
//       orientation of molecule should not dipend on the position of mass of centre
         Rpt[id]  = 0.0;
         RCOM[id] = 0.0;
      }

      Eulang[PHI]=MCAngles[PHI][it_rot];
      Eulang[CTH]=acos(MCAngles[CTH][it_rot]);
      Eulang[CHI]=MCAngles[CHI][it_rot];

      vcord_(Eulang,RCOM,Rpt,vtable,&Rgrd,&THgrd,&CHgrd,&Rvmax,&Rvmin,&Rvstep,&vpot3d,&radret,&theret,&chiret,hatx,haty,hatz,&ivcord);

/*
// get the orientation of the dopant in the next time slice
      int it_rot2 = ((it0 + 1) % NumbTimes)/RotRatio;

      it_rot2 += MCAtom[IMTYPE].offset;

      double hatx2[3];
      double haty2[3];
      double hatz2[3];

      Eulang[PHI]=MCAngles[PHI][it_rot2];
      Eulang[CTH]=acos(MCAngles[CTH][it_rot2]);
      Eulang[CHI]=MCAngles[CHI][it_rot2];

      vcord_(Eulang,RCOM,Rpt,vtable,&vpot3d,&radret,&theret,&chiret,hatx2,haty2,hatz2,&ivcord);
*/

/*
//    check hatz with MCCosine
      for(int id=0;id<NDIM;id++)
      {
         if(fabs(hatz[id] - MCCosine[id][it_rot] > 0.00005))
         cout<<id<<" "<<it_rot<<" "<<hatz[id]<<" "<<MCCosine[id][it_rot]<<endl;
      }
*/
      
//    SFF axes
      if(iframe == 0)
      {
         hatx[0]=1.0,hatx[1]=0.0,hatx[2]=0.0;
         haty[0]=0.0,haty[1]=1.0,haty[2]=0.0;
         hatz[0]=0.0,hatz[1]=0.0,hatz[2]=1.0;
      }

      for (int atom=0;atom<MCAtom[BSTYPE].numb;atom++)
      {

         int it1 = (it0 + 1) % NumbTimes;

         int pt0 = offset + NumbTimes*atom;     // offset only
         int pt1 = pt0;

         if (it1!= (it0 + 1))
         pt1 = offset + NumbTimes*PIndex[atom]; // offset only

         pt0 += it0;
         pt1 += it1;

         for (int dim=0;dim<NDIM;dim++)         // COM adjustment 
         {
            dr0[dim] = MCCoords[dim][pt0] - newcoords[dim][it0];
            dr1[dim] = MCCoords[dim][pt1] - newcoords[dim][it1];
         }

//       get the x0, y0, z0, x1, y1, and z1 in their respective dopant-fixed frame
/*
         double x0=0.0;
         double y0=0.0;
         double z0=0.0;
         double x1=0.0;
         double y1=0.0;
         double z1=0.0;
         for (int dim=0;dim<NDIM;dim++)
         {
            x0 += hatx[dim]*dr0[dim];
            y0 += haty[dim]*dr0[dim];
            z0 += hatz[dim]*dr0[dim];
            x1 += hatx2[dim]*dr1[dim];
            y1 += haty2[dim]*dr1[dim];
            z1 += hatz2[dim]*dr1[dim];
         }

         area_proj[AXIS_X] += 0.5*(y0*z1 - z0*y1);
         area_proj[AXIS_Y] += 0.5*(z0*x1 - x0*z1);
         area_proj[AXIS_Z] += 0.5*(x0*y1 - y0*x1);
*/
//----------- cross product of the two adjacent position vectors with respect to the dopant COM to tet area --------------------------------------

///*
         area[AXIS_X]  = 0.5*(dr0[AXIS_Y]*dr1[AXIS_Z] - dr0[AXIS_Z]*dr1[AXIS_Y]);
         area[AXIS_Y]  = 0.5*(dr0[AXIS_Z]*dr1[AXIS_X] - dr0[AXIS_X]*dr1[AXIS_Z]);
         area[AXIS_Z]  = 0.5*(dr0[AXIS_X]*dr1[AXIS_Y] - dr0[AXIS_Y]*dr1[AXIS_X]);

//       project the vector area onto the three principal axes of the dopant
         for (int id=0;id<NDIM;id++)
         {
            area_projx += area[id] * hatx[id];
            area_projy += area[id] * haty[id];
            area_projz += area[id] * hatz[id];
         }
//*/

///*     calculate the classical moment of inertia
         for (int id=0;id<NDIM;id++)
         {
//          diagonal elements first
            double *hat_dum;
            if(id == AXIS_X)
            hat_dum = hatx;

            if(id == AXIS_Y)
            hat_dum = haty;

            if(id == AXIS_Z)
            hat_dum = hatz;

//          cross product dr0 x hat(x,y,z)
            rn0[AXIS_X] = dr0[AXIS_Y]*hat_dum[AXIS_Z] - dr0[AXIS_Z]*hat_dum[AXIS_Y];
            rn0[AXIS_Y] = dr0[AXIS_Z]*hat_dum[AXIS_X] - dr0[AXIS_X]*hat_dum[AXIS_Z];
            rn0[AXIS_Z] = dr0[AXIS_X]*hat_dum[AXIS_Y] - dr0[AXIS_Y]*hat_dum[AXIS_X];

//          cross product dr1 x hat(x,y,z)
            rn1[AXIS_X] = dr1[AXIS_Y]*hat_dum[AXIS_Z] - dr1[AXIS_Z]*hat_dum[AXIS_Y];
            rn1[AXIS_Y] = dr1[AXIS_Z]*hat_dum[AXIS_X] - dr1[AXIS_X]*hat_dum[AXIS_Z];
            rn1[AXIS_Z] = dr1[AXIS_X]*hat_dum[AXIS_Y] - dr1[AXIS_Y]*hat_dum[AXIS_X];

//          accumulate the dot product of rn0 and rn1 into the diagonal element of classical moment of inertia
            double sum=0.0;
            for (int dim=0;dim<NDIM;dim++)
            sum += rn0[dim]*rn1[dim]*bmass;
//          inert3D[id*NDIM+id] += rn0[dim]*rn1[dim]*bmass;

            if(id == AXIS_X)
            icxx += sum;

            if(id == AXIS_Y)
            icyy += sum;

            if(id == AXIS_Z)
            iczz += sum;

//          off diagonal elements below
//          projection of dr0 on hat_dum (id)
            double dr0_id = 0.0;
            for (int dim=0;dim<NDIM;dim++)
            dr0_id +=hat_dum[dim]*dr0[dim];

            for (int jd=0;jd<NDIM;jd++)
            if( jd != id )
            {
               double *hat_dum2;
               if(jd == AXIS_X)
               hat_dum2 = hatx;

               if(jd == AXIS_Y)
               hat_dum2 = haty;

               if(jd == AXIS_Z)
               hat_dum2 = hatz;

               double dr1_jd = 0.0;
               for (int dim=0;dim<NDIM;dim++)
               dr1_jd +=hat_dum2[dim]*dr1[dim];

               //inert3D[id*NDIM + jd] += -bmass * dr0_id * dr1_jd;

               if(id == AXIS_X && jd == AXIS_Y)
               icxy += -bmass * dr0_id * dr1_jd;

               if(id == AXIS_X && jd == AXIS_Z)
               icxz += -bmass * dr0_id * dr1_jd;

               if(id == AXIS_Y && jd == AXIS_X)
               icyx += -bmass * dr0_id * dr1_jd;

               if(id == AXIS_Y && jd == AXIS_Z)
               icyz += -bmass * dr0_id * dr1_jd;

               if(id == AXIS_Z && jd == AXIS_X)
               iczx += -bmass * dr0_id * dr1_jd;

               if(id == AXIS_Z && jd == AXIS_Y)
               iczy += -bmass * dr0_id * dr1_jd;

            }

         }
//*/

      }
   }

// transfer the parallelly accumulated results in the corresponding arrays
   area_proj[AXIS_X] = area_projx;
   area_proj[AXIS_Y] = area_projy;
   area_proj[AXIS_Z] = area_projz;
   inert3D[AXIS_X*NDIM+AXIS_X] = icxx;
   inert3D[AXIS_X*NDIM+AXIS_Y] = icxy;
   inert3D[AXIS_X*NDIM+AXIS_Z] = icxz;
   inert3D[AXIS_Y*NDIM+AXIS_X] = icyx;
   inert3D[AXIS_Y*NDIM+AXIS_Y] = icyy;
   inert3D[AXIS_Y*NDIM+AXIS_Z] = icyz;
   inert3D[AXIS_Z*NDIM+AXIS_X] = iczx;
   inert3D[AXIS_Z*NDIM+AXIS_Y] = iczy;
   inert3D[AXIS_Z*NDIM+AXIS_Z] = iczz;

// scale the inert3D by 1/NumbTimes and add it to the block accumulation
   for (int id=0;id<NDIM*NDIM;id++)
   {
      if(iframe == 1)
      _inert3DMFF[id] += inert3D[id]/(double)NumbTimes;

      if(iframe == 0)
      _inert3DSFF[id] += inert3D[id]/(double)NumbTimes;

   }

// add area_proj to the block accumulation
   int ind = 0;
   for (int id=0;id<NDIM;id++)
   {
//    cout<<area_proj[id]<<BLANK;
      for (int jd=0;jd<=id;jd++)
      {
         if(iframe == 1)
         _areas3DMFF[ind] += area_proj[id]*area_proj[jd];

         if(iframe == 0)
         _areas3DSFF[ind] += area_proj[id]*area_proj[jd];


         ind ++;
      }
   }
// cout<<endl;

}

void SaveAreaEstimators (const char fname [], double acount, long int blocknumb)
{
  const char *_proc_=__func__;    //  SaveAreaDensities() 
 
//------  open file ---------------

  fstream fid;
  string  fsuper;

  fsuper  = fname;
  fsuper += IO_EXT_SUP;

  fid.open(fsuper.c_str(),ios::app | ios::out); io_setout(fid);

  if (!fid.is_open())
 _io_error(_proc_,IO_ERR_FOPEN,fsuper.c_str());

//--------------------------------
//  shift the center of mass 
//--------------------------------

  double mass   = MCAtom[BSTYPE].mass; 
  double lambda = MCAtom[BSTYPE].lambda; 

  double norm   = 2.0/(MCBeta*lambda);

  fid << setw(IO_WIDTH_BLOCK) << blocknumb << BLANK;     // block number

// super density (normalization: inertia does not include mass)

  fid << setw(IO_WIDTH) << _area2[PERP]*norm/_inert[PERP] << BLANK; 
  fid << setw(IO_WIDTH) << _area2[PARL]*norm/_inert[PARL] << BLANK;

// moment of inertia

  double units = Units.mass*Units.length*Units.length;

  fid << setw(IO_WIDTH) << _inert[PERP]*units*mass/acount << BLANK;
  fid << setw(IO_WIDTH) << _inert[PARL]*units*mass/acount << BLANK;

// area (normalized and averaged)

  fid << setw(IO_WIDTH) << _areas[PERP]*sqrt(norm/_inert[PERP])/acount << BLANK; 
  fid << setw(IO_WIDTH) << _areas[PARL]*sqrt(norm/_inert[PARL])/acount << BLANK;
 
  fid << endl;
  fid.close();
}

void super_reset (void)
{
   _areas[PERP] = 0.0;
   _areas[PARL] = 0.0;

   _area2[PERP] = 0.0;
   _area2[PARL] = 0.0;

   _inert[PERP] = 0.0;
   _inert[PARL] = 0.0;

   for (int id=0;id<6;id++)
   {
      _areas3DMFF[id] = 0.0;
      _areas3DSFF[id] = 0.0;
   }

   for (int id=0;id<NDIM*NDIM;id++)
   {
      _inert3DMFF[id]=0.0;
      _inert3DSFF[id]=0.0;
   }

}

void SaveAreaEstim3D (const char fname [], double acount, long int blocknumb, int iframe)
{

  const char *_proc_=__func__;    //  SaveAreaDensities()

//------  open file ---------------

  fstream fid;
  string  fsuper;

  fsuper  = fname;

  if(iframe == 1)
  fsuper += IO_EXT_MFFSUP3D;

  if(iframe == 0)
  fsuper += IO_EXT_SFFSUP3D;

  fid.open(fsuper.c_str(),ios::app | ios::out); io_setout(fid);

  if (!fid.is_open())
 _io_error(_proc_,IO_ERR_FOPEN,fsuper.c_str());

  double bmass = MCAtom[BSTYPE].mass;
  double lambda = MCAtom[BSTYPE].lambda;

  double norm = 2.0*bmass/(MCBeta * lambda);  // norm = 4m^2/(hbar^2 * beta)

  fid << setw(IO_WIDTH_BLOCK) << blocknumb << BLANK;     // block number

// 9 components of classical moment of inertia
  for(int id=0;id<NDIM*NDIM;id++)
  {
     if(iframe == 1)
     fid << setw(IO_WIDTH) << _inert3DMFF[id]/(double)acount<< BLANK;

     if(iframe == 0)
     fid << setw(IO_WIDTH) << _inert3DSFF[id]/(double)acount<< BLANK;
  }

// 6 components of 4m^2/(hbar^2*beta) <A_iA_j>
   for(int id=0;id<6;id++)
   {
      if(iframe == 1)
      fid << setw(IO_WIDTH) << _areas3DMFF[id]*norm/(double)acount<< BLANK;

      if(iframe == 0)
      fid << setw(IO_WIDTH) << _areas3DSFF[id]*norm/(double)acount<< BLANK;
   }

// 3 compponents of area vector
// for(int id=0;id<NDIM;id++)
// fid << setw(IO_WIDTH) << _area3D[id]/(double)acount <<BLANK;

// fid << setw(IO_WIDTH) << norm <<BLANK;

  fid << endl;
  fid.close();

}

void GetPermutation()
{

   const char *_proc_=__func__;    //  SaveAreaDensities()

   fstream fid;
   string fname;
   fname = FPERMU;

   fid.open(fname.c_str(),ios::app | ios::out); io_setout(fid);

   if (!fid.is_open())
   _io_error(_proc_,IO_ERR_FOPEN,fname.c_str());

   for(int iat=0;iat<MCAtom[BSTYPE].numb;iat++)
   {
      fid<<PIndex[iat]<<" ";
   }

   fid<<endl;
   fid.close();

}
