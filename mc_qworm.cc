#include "mc_qworm.h"
#include "mc_setup.h"
#include "mc_utils.h"
#include "mc_confg.h"
#include "mc_randg.h"
#include "mc_piqmc.h"

#include <stdlib.h>
#include <math.h>

// TEST
const int    MAXNEIGHBORS = 100;    // max number of atoms in the list of neighbors
const double MC_CUTOFF    = 100.0;  // cuttoff for the list of neigbors in thermal wavelength  

TPathWorm Worm;

void  qworm_close  (void);
void  qworm_open   (void);
void  qworm_remove (void);
void  qworm_insert (void);
void  qworm_advance(void); 
void  qworm_recede (void); 
void  qworm_swap   (void);

double qw_open_prob(int);

void   sample_middle (int,int,int,int,double **);

double get_potential (int,int,int,int,double **);

double **QWTotal;    // total number of the worm moves
double **QWAccep;    // total number of the accepted worm moves

double countQW;      // total count of the worm moves

//--- IMPLEMENTATION ------------------------------------

int get_ptable(int,int,int,int,int);
int atom2swap(int, double &);  

double *  _dr2_list;
int    *  _atm_list;

double * _ptable;
double mc_cutoff2;

double _qw_norm;    // normalization for open/close moves
 
void MCWormInit(void)
{
   const char *_proc_=__func__;    // "MCWormInit()";
 
   Worm.exists = 0;   // start with Z-space 

// get a particle type (int) for the worm algorithm

   string stype = Worm.stype;
   bool found   = false;

   for (int type=0;type<NumbTypes;type++)
   if  (stype == MCAtom[type].type)
   {
       Worm.type = type;
       found     = true;  
       break;
   } 
 
   if (!found)
   nrerror(_proc_,"Can't find a particle type for the worm algorithm");

   int type = Worm.type;

   Worm.c *= (Density/(MCAtom[type].numb*NumbTimes*Worm.m));
  _qw_norm = Worm.c  * MCAtom[type].numb*NumbTimes*Worm.m;

  _dr2_list =  new double [MCAtom[type].numb];
  _atm_list =  new int    [MCAtom[type].numb];
  _ptable   =  new double [MCAtom[type].numb+1];

   double termal2 = (double)Worm.m*4.0*MCAtom[type].lambda*MCTau;
   mc_cutoff2 = MC_CUTOFF*MC_CUTOFF*termal2; // the thermal wave length^2 for the segment
}

void MCWormDone(void)
{
  const char *_proc_=__func__;  // "MCWormDone()";

  delete [] _dr2_list;
  delete [] _atm_list;
  delete [] _ptable;
}

void MCWormMove(void)
{
// try many moves before evaluating the average
// for (int time=0;time<NumbTimes;time++)             // loop over segments' origins
for (int atom=0;atom<MCAtom[Worm.type].numb;atom++)   // one atom to move only
{
    countQW += 1.0;

    if (Worm.exists)
    qworm_close();
    else
    qworm_open ();

    if (Worm.exists)
    {
       countQW += 1.0;

       double r = rnd5();
     
       if (r>0.5) 
       qworm_advance(); 
       else
       qworm_recede(); 
    }

    if ((BOSONS) && (Worm.type == BSTYPE))
    {
      countQW += 1.0;
      if (Worm.exists) 
      qworm_swap();
    }
} // END loop over atom labels 
}

double qw_open_prob(int segm)
{
   int offset = MCAtom[Worm.type].offset;

   double kin = 0.0;          // free particle density (m*\tau) 

   int pt0 = offset + Worm.atom_i*NumbTimes + Worm.ira;        
   int pt1 = offset + Worm.atom_m*NumbTimes + Worm.masha;

   for (int id=0;id<NDIM;id++)
   {
      double dr = MCCoords[id][pt0] - MCCoords[id][pt1];

      if (MINIMAGE)
      dr  -= (BoxSize*rint(dr/BoxSize));

      kin += (dr*dr); 
   }

   kin /= (MCAtom[Worm.type].twave2*(double)segm);    // 4\lambda m \tau

   double pot = get_potential(Worm.ira,Worm.ira+segm,Worm.atom_i,Worm.atom_m,MCCoords);
   pot *= MCTau;

// return (_qw_norm*sqrt((double)segm)*exp(kin+pot));
   return (_qw_norm*pow((double)segm,0.5*(double)NDIM)*exp(kin+pot));
}

void  qworm_open(void)
{
   QWTotal[0][QW_OPEN] += 1.0;

   Worm.atom_i = nrnd1(MCAtom[Worm.type].numb); // select a world line
   Worm.ira    = nrnd2(NumbTimes);              // select a time slice
   int segm    = nrnd3(Worm.m) + 1;             // select number of slices to remove

   Worm.masha  = (Worm.ira + segm) % NumbTimes;
   Worm.atom_m = Worm.atom_i;

   if (Worm.masha != (Worm.ira + segm))        // do not need this if and the previous line
   Worm.atom_m = PIndex[Worm.atom_i];          // if PIndex is defined for Boltzmann statistics

   double prob = qw_open_prob(segm);

   bool Accepted = false;
   if (prob >=1.0)          Accepted = true;
   else if (prob > rnd1())  Accepted = true;

   if (Accepted)     
   {
      Worm.exists          = 1;
      QWAccep[0][QW_OPEN] += 1.0;
   } 
}

void  qworm_close  (void)
{
   QWTotal[0][QW_CLOSE] += 1.0;

   int segm  = Worm.masha - Worm.ira;

   if (segm<0)
   segm  += NumbTimes;

   if (segm>Worm.m) return;  // should be balanced with the "open" move

#ifdef DEBUG_WORM
   if (segm == 0)  nrerror("qworm_close: ", "ira - masha = 0"); 
#endif

   int it0 = Worm.ira;
   int it2 = Worm.ira + segm; // masha = it2 % NumbTimes

/* 
   int offset = MCAtom[Worm.type].offset;

   int pt0 = offset + Worm.atom_i*NumbTimes + Worm.ira;        
   int pt1 = offset + Worm.atom_m*NumbTimes + Worm.masha;

// kinetic 

   double kin = 0.0;          // free particle density (m*\tau) 
   for (int id=0;id<NDIM;id++)
   {
      double dr = MCCoords[id][pt0] - MCCoords[id][pt1];
 
      if (MINIMAGE)
      dr  -= (BoxSize*rint(dr/BoxSize));

      kin += (dr*dr); 
   }

   kin /= (MCAtom[Worm.type].twave2*(double)segm);    // 4\lambda m \tau

   if (exp(-kin) < rnd1()) return;                    // normalization ? 
*/ 
   sample_middle(it0,it2,Worm.atom_i,Worm.atom_m,MCCoords);

   double prob = 1.0/qw_open_prob(segm);

   bool Accepted = false;
   if (prob >=1.0)          Accepted = true;
   else if (prob > rnd1())  Accepted = true;
  
   if (Accepted)   // remove a worm 
   {
      Worm.exists = 0; 
      QWAccep[0][QW_CLOSE] += 1.0;
   }
}

void  sample_middle(int it0, int it2, int atom0, int atom2, double ** coords)
//
//   it2-it0 should be short enough (< NumbTimes) :  p[atom0] = atom2 
//                                     
{
   if ((it2-it0)<2) return;

   int it1 = (int)rint(0.5*(double)(it0+it2)); // the middle point

   int pt0 = it0 % NumbTimes;
   int pt1 = it1 % NumbTimes;
   int pt2 = it2 % NumbTimes;

   int atom1 = atom0;

   if ((pt1 != it1) && (pt0 == it0)) 
   {
       atom1 = atom2;   
#ifdef DEBUG_WORM
       if ((atom1 != PIndex[atom0]) || (atom0 != RIndex[atom1]))
       nrerror("sample_middle(): ","Wrong permutation type"); 
#endif
   }

   int offset = MCAtom[Worm.type].offset;

   pt0 += (offset + atom0*NumbTimes);
   pt1 += (offset + atom1*NumbTimes);
   pt2 += (offset + atom2*NumbTimes);

   double s0 = (double)(it1-it0);  
   double s2 = (double)(it2-it1);  

   double gkin = (s0+s2)/(MCAtom[Worm.type].twave2*s0*s2); 

   for (int id=0;id<NDIM;id++)
   {  	  
      coords[id][pt1]  = (s2*coords[id][pt0]+s0*coords[id][pt2])/(s0+s2);
      coords[id][pt1] += gauss(gkin);
   } 

   sample_middle(it0,it1,atom0,atom1,coords);
   sample_middle(it1,it2,atom1,atom2,coords);

   return;
}

void  qworm_remove (void)
{
   return;
}

void  qworm_insert (void)
{
   return;
}

void  qworm_advance(void)
{
// const char *_proc_ = __func__;  // qworm_advance()

   QWTotal[0][QW_ADVANCE] += 1.0;

#ifdef DEBUG_WORM 
   if (Worm.ira == Worm.masha)    
   nrerror("qworm_advance: ", "Worm.ira - Worm.masha = 0"); 
#endif
 
   int segm  = Worm.masha - Worm.ira;

   if (segm<0)
   segm  += NumbTimes;

   int advance  = nrnd3(Worm.m) + 1; //#beads = #tslices

   if (segm-advance<=0) return; // canonical, do not introduce a new particle

   int type   = Worm.type;
   int offset = MCAtom[type].offset;

// kinetic contribution : sample a segment
   
   int it0 = Worm.ira;
   int it2 = Worm.ira + advance;
  
   int ira_new     = it2 % NumbTimes;
   int atom_i_new  = Worm.atom_i;

   if (ira_new != it2)
   atom_i_new = Worm.atom_m;

// sample a new position for the end point

   double gvar = 1.0/((double)advance*MCAtom[type].twave2); // variance for gaussian sampling 

   int pt0 = offset + Worm.atom_i*NumbTimes + it0 % NumbTimes;
   int pt2 = offset + atom_i_new *NumbTimes + it2 % NumbTimes;

   for (int id=0;id<NDIM;id++)
   MCCoords[id][pt2] = MCCoords[id][pt0] + gauss(gvar);

   sample_middle(it0,it2,Worm.atom_i,atom_i_new,MCCoords);

   double pot = get_potential(it0,it2+1,Worm.atom_i,atom_i_new,MCCoords); // no 0.5 Pot !

   bool Accepted = false;
   if (pot < 0.0)             Accepted = true;
   else if
   (exp(-pot*MCTau) > rnd2()) Accepted = true;

   if (Accepted)           
   {
      QWAccep[0][QW_ADVANCE] += 1.0;

      Worm.ira    = ira_new;
      Worm.atom_i = atom_i_new; // Worm.atom_m
   }
} 

void  qworm_recede (void)
{
// const char *_proc_ = __func__;  // qworm_recede()

   QWTotal[0][QW_RECEDE] += 1.0;
 
   int segm   = Worm.ira - Worm.masha;

   if (segm < 0)
   segm  += NumbTimes;

   int recede = nrnd3(Worm.m) + 1; 

   if ((segm - recede) < 1) return;  // canonical, do not remove a particle

   int it0   = (Worm.ira - recede);
   int it1   =  Worm.ira;

   int atom0 = Worm.atom_i;
   int atom1 = Worm.atom_i;

   if (it0<0)
   {
      it0 += NumbTimes;
      it1 += NumbTimes; 
 
      atom0 = RIndex[atom1];   // canonical simulations
   }

   double pot = get_potential(it0,it1+1,atom0,atom1,MCCoords); // no 0.5 contribution in Pot

   bool Accepted = false;
   if (pot > 0.0)           Accepted = true;
   else if
  (exp(pot*MCTau) > rnd2()) Accepted = true;

   if (Accepted)           
   {
      Worm.ira    = it0 % NumbTimes;
      Worm.atom_i = atom0;
      
      QWAccep[0][QW_RECEDE] += 1.0; 
   }
}
 
double get_potential(int it0, int it1, int atom0, int atom1, double ** coords)
{
   int atom_offset = MCAtom[Worm.type].offset/NumbTimes; // global numeration of atoms in PotEnergy() 
   
   int pit0 = it0 % NumbTimes;  // atom0 it's important for some functions 
   int pit1 = it1 % NumbTimes;  // atom1

   double pot = 0.0;

   int atom = atom0;
   for (int it=(it0+1);it<it1;it++)
   {
       int pit = it % NumbTimes;

       if ((pit != it) && (pit0 == it0)) // could happen only once
       atom = atom1;                     // should be equivalent to PIndex[atom0]

       pot += PotEnergy(atom_offset + atom, coords,pit);   
   }
   return pot;
}

void  qworm_swap   (void)
{
   QWTotal[0][QW_SWAP] += 1.0;

   int segm  = Worm.m;

   int  it0  = Worm.ira;     
   int  it1  = it0 + segm; 

   int pit0  = it0;    // [r(pt0)-r(pt1)]^2
   int pit1  = it1 % NumbTimes; 

   int atomw = Worm.atom_i;

   int count = get_ptable(atomw,pit0,pit1,segm,it1); // number of atoms in the neighborhood 

   if (count<=0) return;                             // empty permutation table

   double pnorm_old;
   double pnorm_new;
 
   int atom1 = atom2swap(count,pnorm_old);           // guess for swaping 
  
   if (atom1 < 0) return; // identity permuation, do nothing
 
   int atom0 = atom1;
   if (pit1 != it1)
   atom0 = RIndex[atom1];

#ifdef DEBUG_WORM
   if (atom0 == Worm.atom_i)   // masha between it0 and it1
   nrerror("qworm_swap: ","Wrong permutation type 1");   

   if ((Worm.masha < Worm.ira) && (Worm.atom_i == Worm.atom_m) && (atom1 == Worm.atom_i))
   nrerror("qworm_swap: ","Wrong permutation type 2");   
#endif

   int type    =  Worm.type;

   int offset0 = MCAtom[type].offset + atom0*NumbTimes;
   int offset1 = MCAtom[type].offset + atom1*NumbTimes;
   int offsetw = MCAtom[type].offset + atomw*NumbTimes;

   for (int id=0;id<NDIM;id++) // init the end points
   { 
      newcoords[id][offset0 + pit0] = MCCoords[id][offsetw + pit0]; 
      newcoords[id][offset1 + pit1] = MCCoords[id][offset1 + pit1]; 
   }

   sample_middle(it0,it1,atom0,atom1,newcoords);

   int gatom0 = offset0/NumbTimes; // offset/NumbTimes 
   int gatom1 = offset1/NumbTimes; // offset/NumbTimes 

   double pot = 0.0;
   int gatom = gatom0;

   for (int it=(it0+1);it<it1;it++)
   {
      int pit = it % NumbTimes;

      if (pit != it)
      gatom = gatom1;
    
      pot += PotEnergy(gatom, newcoords, pit); // new configs
      pot -= PotEnergy(gatom, MCCoords,  pit); // old configs
   }

   double prob = exp(-pot*MCTau);

   count = get_ptable(atom0,pit0,pit1,segm,it1); 
   pnorm_new = 0.0;
   for (int ic=1;ic<=count;ic++) pnorm_new += _ptable[ic];

   prob *= (pnorm_old/pnorm_new); 

   bool Accepted = false;
   if (prob >=1.0)          Accepted = true;
   else if (prob > rnd4())  Accepted = true;
 
   if (Accepted)           
   {
      QWAccep[0][QW_SWAP] += 1.0; 

//    int offset = offset0;      

      for (int id=0;id<NDIM;id++)
      {
 
      int offset = offset0;      

      for (int it=(it0+1);it<it1;it++) 
      {
         int pit = it % NumbTimes;
         if (pit != it)
         offset = offset1;

         MCCoords[id][offset + pit] = newcoords[id][offset + pit]; 
      }
      } 

      for (int id=0;id<NDIM;id++)  
      for (int it=0;it<=it0;it++)   
      {
         newcoords[id][offset0 + it] = MCCoords [id][offset0 + it]; // save
         MCCoords [id][offset0 + it] = MCCoords [id][offsetw + it]; // swap 
         MCCoords [id][offsetw + it] = newcoords[id][offset0 + it]; // swap
      }

// SWAP: permutation table  ----------------

      int ratomw = RIndex[atomw];  
      int ratom0 = RIndex[atom0];

      PIndex[ratomw] =  atom0;
      RIndex[atom0 ] =  ratomw;

      PIndex[ratom0] =  atomw;
      RIndex[atomw]  =  ratom0;

      if (Worm.ira > Worm.masha)
      if (atom0 == Worm.atom_m)
      Worm.atom_m = Worm.atom_i; 
      else   
      if (Worm.atom_i == Worm.atom_m)
      Worm.atom_m = atom0;
   }  // end Accepted
} 

bool WorldLine(int atom, int pt)
//
//  true   if pt belongs to the world line
//  false  if pt belongs to the gap between masha and ira
//
{
    bool wline = true; 

    if ((atom == Worm.atom_m) || (atom == Worm.atom_i))
    if ((Worm.atom_i !=  Worm.atom_m)||(Worm.ira > Worm.masha)) // Worm.ira > Worm.masha or swap test 
    { 
       if (((atom == Worm.atom_m) && (pt < Worm.masha))
        || ((atom == Worm.atom_i) && (pt > Worm.ira)))
        wline = false;
    } 
    else // Worm.ira < Worm.masha)    // DO NOT NEED IF Worm.atom_i !=  Worm.atom_m 
    {
        if ((pt > Worm.ira) && (pt < Worm.masha))
        wline = false;
    }

    return wline;
}

int get_ptable(int atomw, int pt0, int pt1, int segm, int t1)
// RETURN: number of entries in the permutation table
//         ptable [1...count]
{
  int type   = Worm.type;
  int offset = MCAtom[type].offset;
  int itw    = offset + atomw*NumbTimes + pt0;

  int count = 0;   // number of atoms in the list of neighbours

  for (int atom1=0;atom1<MCAtom[type].numb;atom1++) 
  if (WorldLine(atom1,pt1))
  {
     int atom0  = atom1;
    
     if  (t1 != pt1)
     {
        atom0 = RIndex[atom1];
#ifdef DEBUG_WORM
        if (PIndex[atom0] != atom1)   // DEBUG
        {
           cout <<"ERROR IN PT"<<endl;
           exit(0);   
        } 
#endif
     }
 
     if( atom0 != Worm.atom_i)  
     { 
        int it1 = offset + atom1*NumbTimes + pt1;
        
        double dr2 = 0.0; 

        for (int id=0;id<NDIM;id++)
        {
           double dx = MCCoords[id][itw] - MCCoords[id][it1];  // exchange

           if (MINIMAGE)
           dx  -= (BoxSize*rint(dx/BoxSize));

           dr2 += (dx*dx);
        }

        if (dr2 < mc_cutoff2) // insert an atom in the list of neighbours
        {
           count++;        
          _dr2_list [count] = dr2;    // squared distances
          _atm_list[count]  = atom1;
        }
    }
  } // end the loop over atom1
  
  mmsort(_dr2_list,_atm_list,count);  // count:  [1, count]

  if (count>MAXNEIGHBORS)  //cutoff in the permutation space 
  count = MAXNEIGHBORS;
      
// build the permutation table

  double norm = 1.0/((double)segm*MCAtom[type].twave2);

  for (int ic=1;ic<=count;ic++)
 _ptable[ic] = exp(-norm*_dr2_list[ic]);

  return count;
}

int atom2swap(int count, double & pnorm)  // pick the atom to swap
{
// treat count == 1 case separately
// if (count <= 0)  // check outside       
// nrerror("atom2swap","empty list of neighbors");

   pnorm = 0.0;
   for (int ic=1;ic<=count;ic++) pnorm += _ptable[ic];
 
   double prand = pnorm*rnd3();

   double sum = 0.0;
   int ic     = 1;
   while ((ic<=count) && (sum<prand)) // pick the permutation
   {
      sum += _ptable[ic];
      ic++;
   }  
//------------------------- it's important
   ic--;   
//-------------------------
  
   return (_atm_list[ic]); // return the atom label to swap 
}

void ResetQWCounts(void)
{
   countQW = 1.0; 

   for (int im=0;im<QWMAXMOVES;im++)
   {
      QWTotal[0][im] = 0.0;
      QWAccep[0][im] = 0.0;
   }
}

void MemAllocQWCounts(void)
{
   QWTotal = doubleMatrix(1,QWMAXMOVES);
   QWAccep = doubleMatrix(1,QWMAXMOVES);
}

void MFreeQWCounts(void)
{
   free_doubleMatrix(QWTotal);
   free_doubleMatrix(QWAccep);
}
