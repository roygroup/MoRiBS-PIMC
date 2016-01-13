#include "mc_poten.h"
#include "mc_confg.h"
#include "mc_input.h"
#include "mc_utils.h"
#include "mc_setup.h"

#include <math.h>
#include <stdlib.h>
#include <string>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <stdio.h>

// ----------- POTENTIALS ----------------------------------
// 1D
double * _poten1D [MAX_NUMBER_INTER]; // tabulated 1D potentials
double * _pgrid1D [MAX_NUMBER_INTER]; // grids for 1D potentials
double * _pderiv2 [MAX_NUMBER_INTER]; // second derives for spline

int      _psize1D [MAX_NUMBER_INTER]; // number of grid points for 1D potential

double   _alpha[MAX_NUMBER_INTER];    // parameters for u0e^(-alpha*r) extrapolation
double   _unode[MAX_NUMBER_INTER];    // of the potential
double   _c6   [MAX_NUMBER_INTER];    // -C6/r^6

// 2D revised by Hui Li
//double ** _dpoten2D [MAX_NUMBER_INTER]; // tabulated 2D potentials difference between ground and excited states (V1-V0)
double *  _drgrid2D [MAX_NUMBER_INTER]; // r    grid points for 2D
double *  _dcgrid2D [MAX_NUMBER_INTER]; // cost grid points for 2D

int       _drsize2D [MAX_NUMBER_INTER]; // number of r   grid points  for 2D potential
int       _dcsize2D [MAX_NUMBER_INTER]; // number of cost grid points for 2D potential

double    _ddelta_r [MAX_NUMBER_INTER]; //  delta r    for 2D potential
double    _ddelta_c [MAX_NUMBER_INTER]; //  delat cost for 2D potential


//2D 
double ** _poten2D [MAX_NUMBER_INTER]; // tabulated 2D potentials
double *  _rgrid2D [MAX_NUMBER_INTER]; // r    grid points for 2D
double *  _cgrid2D [MAX_NUMBER_INTER]; // cost grid points for 2D

int       _rsize2D [MAX_NUMBER_INTER]; // number of r   grid points  for 2D potential
int       _csize2D [MAX_NUMBER_INTER]; // number of cost grid points for 2D potential

double    _delta_r [MAX_NUMBER_INTER]; //  delta r    for 2D potential
double    _delta_c [MAX_NUMBER_INTER]; //  delat cost for 2D potential

// Non-linear rotor, 3D potential
double * vtable ;//= new double [SizePotTab];            //  3D potential table added by Toby
int      Rgrd;  // # of radial grid points to store 3D potential
int      THgrd; // # of theta grid points to store 3D potential. It should be 181 at this moment.
int      CHgrd; // # of chi grid points to store 3D potential. It should be 91, 181, or 361 at this moment.
double   Rvmax; // maximum radius for 3D potential extrapolation
double   Rvmin; // minimum radius for 3D potential extrapolation
double   Rvstep; // radial increment for 3D potential extrapolation

// -------  ROTATIONAL DENSITY MATRICES -------------------------

double * _rotgrid [MAX_NUMBER_ROTDN];    // grid points for the rotational density matrix
double * _rotdens [MAX_NUMBER_ROTDN];    // rotational density matrix
double * _rotderv [MAX_NUMBER_ROTDN];    // derivatives of the rotational density matrix
double * _rotesqr [MAX_NUMBER_ROTDN];    // second derivatives of the rotational density matrix

int      _rotsize [MAX_NUMBER_ROTDN];   // number of grid points for the rot density

double * _rotdens_drv2 [MAX_NUMBER_ROTDN]; // second derivatives for spline
double * _rotderv_drv2 [MAX_NUMBER_ROTDN]; // second derivatives for spline
double * _rotesqr_drv2 [MAX_NUMBER_ROTDN]; // second derivatives for spline

extern "C" void rotred_(double *rhoprp,double *erotpr);

extern "C" void potred_(char *fname,double *vtable);

void init_rotdens(int type);
void init_rot3D(int type);

// --------------------------------------------------------------

int  get_filesize (const char []);
void read_datafile(const char [],double *,double *);
void read_datafile(const char [],double *,double *,double *);
void read_datafile(const char [],double *,double *,double *,double *);

void init_pot1D(int atype);
void init_pot2D(int atype);  
void init_pot3D(int atype);  
//void init_dpot2D(int atype);  //revised by Hui Li

//-----------------------------------------------------------

void InitPotentials(void)
{
   const char *_proc_=__func__;    // "InitPotentials"; 

   for (int atype=0;atype<NumbTypes;atype++)
   if ((MCAtom[atype].molecule == 1)||(MCAtom[atype].molecule == 2))                            // molecules
   {
      if ((MCAtom[atype].molecule == 1) && ((MCAtom[atype].numb  > 1) || (MCAtom[atype].numb  < 0)) )
      nrerror(_proc_,"No more than one linear dopant molecule so far"); // check potential energy

      if ((MCAtom[atype].molecule == 2) && ((MCAtom[atype].numb  > NumbRotLim) || (MCAtom[atype].numb  < 0)) )
      nrerror(_proc_,"Weird # of non-linear rotors"); // check # of non-linear rotors

      if(MCAtom[atype].molecule == 1)
      init_pot2D(atype);
//      init_dpot2D(atype);   //revised by Hui Li

      if(MCAtom[atype].molecule == 2 && NumbTypes > 1)
      {
      init_pot3D(atype);
//    potred_(vtable);
      }

   }
   else                                                  // atoms
   init_pot1D(atype);
}

void DonePotentials(void)
{

   for (int atype=0;atype<NumbTypes;atype++)
   if (MCAtom[atype].molecule == 1)              // molecules
   {
       delete [] _rgrid2D[atype];
       delete [] _cgrid2D[atype];
 
//       delete [] _drgrid2D[atype];     //add by Hui Li  
//       delete [] _dcgrid2D[atype];     //add by Hui Li

    
       free_doubleMatrix(_poten2D[atype]);  // atoms 

//       free_doubleMatrix(_dpoten2D[atype]);  // atoms revised by Hui Li 
   }
   else if (MCAtom[atype].molecule == 2)
   {
      delete [] vtable;
   }
   else
   { 
       delete [] _pgrid1D[atype];
       delete [] _poten1D[atype];
       delete [] _pderiv2[atype];  
   }
}

void InitRotDensity(void)
{
   const char *_proc_=__func__;    // "InitRotDensity()";

// this is only performed for a linear rotor.  The non-linear rotor is treated by hard coding by Toby
   if(MCAtom[IMTYPE].molecule == 1 && RotDenType == 0)
   init_rotdens(IMTYPE);

   if(MCAtom[IMTYPE].molecule == 2 && RotDenType == 0)
   {
      init_rot3D(IMTYPE);
//    rotred_(rhoprp,erotpr);
   }

}

void DoneRotDensity(void)
{  
   for (int atype=0;atype<MAX_NUMBER_ROTDN;atype++)
   {
      delete [] _rotgrid [atype];    
      delete [] _rotdens [atype];    
      delete [] _rotderv [atype];   
      delete [] _rotesqr [atype];   

      delete [] _rotdens_drv2 [atype]; 
      delete [] _rotderv_drv2 [atype]; 
      delete [] _rotesqr_drv2 [atype]; 
   }
}

/*
// add by Hui Li for calculation of vibrational shift

void init_dpot2D(int atype)  // read 2D potential difference from the file of .diff
//
//  POTENTIAL FILE FORMAT: #1  [number r points]  [number cost points]
//                         #2   [delta  r ]        [delta  cost]
//                              [r    grid]     monotonically increase
//                              [cost grid]     monotonically increase 
//                    
//                              dpot2D (r,cost)  matrix at grid points
{
    const char *_proc_=__func__;    // "init_dpot2D"; 

    if (atype>=MAX_NUMBER_INTER)
    nrerror(_proc_,ERR_INDEX_EXCEED);

    string fextn;
    fextn = EXT_DIFF; // ".diff" difference potential  

    string fname    =  (MCAtom[atype].fpot+fextn);   // file name

    ifstream fid(fname.c_str(),ios::in);   

    if (!fid.good())
   _io_error(_proc_,IO_ERR_FOPEN,fname.c_str());

//---- read  grid information -------------------------------------
    int drsize;
    int dcsize;

    fid >> drsize;
    fid >> dcsize;
  
   _drsize2D[atype] = drsize;
   _dcsize2D[atype] = dcsize;
    
    fid >> _ddelta_r[atype];
    fid >> _ddelta_c[atype];

   _drgrid2D[atype] = new double [drsize];
   _dcgrid2D[atype] = new double [dcsize];

   _dpoten2D[atype] = doubleMatrix(drsize,dcsize);  

    for (int ir=0;ir<drsize;ir++)  // load r   grid
    fid >> _drgrid2D[atype][ir];         

    for (int ic=0;ic<dcsize;ic++)  // load cost grid
    fid >> _dcgrid2D[atype][ic];      

    for (int ir=0;ir<drsize;ir++)  // read 2D potential        
    for (int ic=0;ic<dcsize;ic++)
    fid >> _dpoten2D[atype][ir][ic];    

    fid.close();

//convert to dimensionless representation 

    drsize = _drsize2D[atype];
    dcsize = _dcsize2D[atype];

   _ddelta_r[atype] /= Units.length;

    for (int ir=0;ir<drsize;ir++)  
   _drgrid2D[atype][ir] /= Units.length;         

    for (int ir=0;ir<drsize;ir++)         
    for (int ic=0;ic<dcsize;ic++)
   _dpoten2D[atype][ir][ic] /=  Units.energy;  
}
*/

void init_pot3D(int atype) // read 3D potential from the file
{

    const char *_proc_=__func__;    // "init_pot2D";

    int pmod = MCAtom[atype].pmod;

    string fextn;
    if      (pmod == PRIMITIVE) fextn = EXT_PRIM; // primitive approximation
    else if (pmod == EFFECTIVE) fextn = EXT_EFFC; // effective potential
    else
    nrerror(_proc_,"Unknown model of interaction");

    string fname    =  (MCAtom[atype].fpot+fextn);   // file name

    char filnam[100];

//  strcpy(filnam,fname.c_str());

    cout<<"potential in "<<fname<<endl;

//  potred_(filnam,vtable);

    ifstream fid(fname.c_str(),ios::in);

    if (!fid.good())
    _io_error(_proc_,IO_ERR_FOPEN,fname.c_str());

//---- read  grid information -------------------------------------
    fid >> Rgrd;
    fid >> THgrd;
    fid >> CHgrd;
    fid >> Rvmin;
    fid >> Rvmax;
    Rvstep=(Rvmax - Rvmin)/(double)(Rgrd-1);
    cout<<"Rgrd="<<Rgrd<<" THgrd="<<THgrd<<" CHgrd="<<CHgrd<<" Rvmin="<<Rvmin<<" Rvmax="<<Rvmax<<" Rvstep="<<Rvstep<<endl;

    if(CHgrd != 361 && CHgrd != 181 && CHgrd != 91 && THgrd != 181)
    nrerror(_proc_,"Degree wise THETA or CHI grid is violated");

    vtable = new double [Rgrd*THgrd*CHgrd];

    for(int line=0;line<Rgrd*THgrd*CHgrd;line++)
    fid >> vtable[line];

    cout<<vtable[0]<<" "<<vtable[Rgrd*THgrd*CHgrd-1]<<endl;
    
    fid.close();

}

void init_pot2D(int atype)  // read 2D potential from the file
//
//  POTENTIAL FILE FORMAT: #1  [number r points]  [number cost points]
//                         #2   [delta  r ]        [delta  cost]
//                              [r    grid]     monotonically increase
//                              [cost grid]     monotonically increase 
//                    
//                              pot2D (r,cost)  matrix at grid points
{
    const char *_proc_=__func__;    // "init_pot2D"; 

    if (atype>=MAX_NUMBER_INTER)
    nrerror(_proc_,ERR_INDEX_EXCEED);

    int pmod = MCAtom[atype].pmod;

    string fextn;
    if      (pmod == PRIMITIVE) fextn = EXT_PRIM; // primitive approximation
    else if (pmod == EFFECTIVE) fextn = EXT_EFFC; // effective potential
    else 
    nrerror(_proc_,"Unknown model of interaction");

    string fname    =  (MCAtom[atype].fpot+fextn);   // file name

    ifstream fid(fname.c_str(),ios::in);   

    if (!fid.good())
   _io_error(_proc_,IO_ERR_FOPEN,fname.c_str());

//---- read  grid information -------------------------------------
    int rsize;
    int csize;

    fid >> rsize;
    fid >> csize;
  
   _rsize2D[atype] = rsize;
   _csize2D[atype] = csize;
    
    fid >> _delta_r[atype];
    fid >> _delta_c[atype];

   _rgrid2D[atype] = new double [rsize];
   _cgrid2D[atype] = new double [csize];

   _poten2D[atype] = doubleMatrix(rsize,csize);  

    for (int ir=0;ir<rsize;ir++)  // load r   grid
    fid >> _rgrid2D[atype][ir];         

    for (int ic=0;ic<csize;ic++)  // load cost grid
    fid >> _cgrid2D[atype][ic];      

    for (int ir=0;ir<rsize;ir++)  // read 2D potential        
    for (int ic=0;ic<csize;ic++)
    fid >> _poten2D[atype][ir][ic]; 

    fid.close();

//convert to dimensionless representation 

    rsize = _rsize2D[atype];
    csize = _csize2D[atype];

   _delta_r[atype] /= Units.length;

    for (int ir=0;ir<rsize;ir++)  
   _rgrid2D[atype][ir] /= Units.length;         

    for (int ir=0;ir<rsize;ir++)         
    for (int ic=0;ic<csize;ic++)
   _poten2D[atype][ir][ic] /=  Units.energy; 
}

void init_pot1D(int atype)
//
//  read/initialize 1D potential MODEL = PRIMITIVE : primitive approximation   
//                               MODEL = EFFECTIVE : effective potential
{
    const char *_proc_=__func__;    // "init_pot1D"; 

    if (atype>=MAX_NUMBER_INTER)
    nrerror(_proc_,ERR_INDEX_EXCEED);

    int pmod = MCAtom[atype].pmod;

    string fextn;
    if      (pmod == PRIMITIVE) fextn = EXT_PRIM; // primitive approximation
    else if (pmod == EFFECTIVE) fextn = EXT_EFFC; // effective potential
    else 
    nrerror(_proc_,"Unknown model of interaction");

    string fname    =  (MCAtom[atype].fpot + fextn);   // file name
    
    int size        =  get_filesize(fname.c_str());  

   _psize1D[atype]  =  size;
   _pgrid1D[atype]  =  new double [size];  // grid
   _poten1D[atype]  =  new double [size];  // data
   _pderiv2[atype]  =  new double [size];  // second derivatives for spline  

    read_datafile(fname.c_str(),_pgrid1D[atype],_poten1D[atype]);

    for (int id=0;id<size;id++)
    {
      _pgrid1D[atype][id] /= Units.length;
      _poten1D[atype][id] /= Units.energy;
    }

    init_spline(_pgrid1D[atype],_poten1D[atype],_pderiv2[atype],_psize1D[atype]);

    if (size<2) nrerror(_proc_,"Not enough data points for extrapolation");

// -------- extrapolation for short distances ------------

//  model pot=U0*exp(-alpha*r) - small distance extrapolation  

    double fr = _poten1D[atype][0]/_poten1D[atype][1];
    double dr = _pgrid1D[atype][1]-_pgrid1D[atype][0];

    _alpha[atype] =  log(fr)/dr;
    _unode[atype] = _poten1D[atype][0]*exp(_alpha[atype]*_pgrid1D[atype][0]);

// -------- extrapolation for long distances -------------

//  model pot=-C6/r^6  - small distance  asymptotic

    double r0 = pow(_pgrid1D[atype][size-2],6.0);
    double r1 = pow(_pgrid1D[atype][size-1],6.0);

    fr = _poten1D[atype][size-1]-_poten1D[atype][size-2];

   _c6[atype] =  fr/(1.0/r0-1.0/r1);
}

void init_rot3D(int type)
{
//  read rotational density matrix ane energy estimator for IMTYPE type asymmetric top molecule
   string  fname = MCAtom[type].type;

   stringstream time; time << NumbRotTimes;                  // number of time slices 
   stringstream temp; temp << Temperature*Units.temperature; // temperature

   fname += ("_T" + temp.str() + "t" + time.str());

   cout<<fname.c_str()<<endl;

   string fden = fname;
   string feng = fname;
   string fesq = fname;

   fden += EXT_RHO;
   feng += IO_EXT_ENG;
   fesq += EXT_ESQ;

   cout<<fden<<" "<<feng<<" "<<fesq<<endl;

   ifstream fid(fden.c_str(),ios::in);

   const char *_proc_=__func__;    // "init_rot3D"

   if (!fid.good())
   _io_error(_proc_,IO_ERR_FOPEN,fden.c_str());

   for(int entry=0;entry<SizeRotDen;entry++)
   fid >> rhoprp[entry];

   fid.close();

   ifstream fid2(feng.c_str(),ios::in);

   if (!fid2.good())
   _io_error(_proc_,IO_ERR_FOPEN,feng.c_str());

   for(int entry=0;entry<SizeRotDen;entry++)
   fid2 >> erotpr[entry];

   fid2.close();

   ifstream fid3(fesq.c_str(),ios::in);

   if (!fid3.good())
   _io_error(_proc_,IO_ERR_FOPEN,fesq.c_str());

// if (!fid3.good())
// {
//    for(int entry=0;entry<SizeRotDen;entry++)
//    erotsq[entry]=0.0;
// }
// else
// {
      for(int entry=0;entry<SizeRotDen;entry++)
      fid3 >> erotsq[entry];

      fid3.close();
// }

// cout<<rhoprp[0]<<" "<<erotpr[0]<<" "<<erotsq[0]<<" "<<rhoprp[SizeRotDen-1]<<" "<<erotpr[SizeRotDen-1]<<" "<<erotsq[SizeRotDen-1]<<endl;

// cout<<erotpr[3999]<<" "<<erotpr[4000]<<endl;

}

void init_rotdens(int type)
//
//  read/initialize rotational density matrix for IMTYPE type molecule
//                             
{
    const char *_proc_=__func__;    // "init_rotdens"; 

    if (type>=MAX_NUMBER_ROTDN)
    nrerror(_proc_,ERR_INDEX_EXCEED);

    string  fname = MCAtom[type].type;

    stringstream time; time << NumbRotTimes;                  // number of time slices 
    stringstream temp; temp << Temperature*Units.temperature; // temperature

    fname += ("_T" + temp.str() + "t" + time.str()); 
    fname += EXT_ROTD;

//  string fname  = (MCAtom[IMTYPE].rdens + (string) EXT_ROTD);   // file name
 
    int size      = get_filesize(fname.c_str());  

   _rotsize[type] = size;
   _rotgrid[type] = new double [size];  
   _rotdens[type] = new double [size];  
   _rotderv[type] = new double [size];  
   _rotesqr[type] = new double [size];  

   _rotdens_drv2[type] = new double [size];  
   _rotderv_drv2[type] = new double [size]; 
   _rotesqr_drv2[type] = new double [size]; 

//  read_datafile(fname.c_str(),_rotgrid[type],_rotdens[type],_rotderv[type]);
    read_datafile(fname.c_str(),_rotgrid[type],_rotdens[type],_rotderv[type],_rotesqr[type]);

    init_spline(_rotgrid[type],_rotdens[type],_rotdens_drv2[type],_rotsize[type]);
    init_spline(_rotgrid[type],_rotderv[type],_rotderv_drv2[type],_rotsize[type]);
    init_spline(_rotgrid[type],_rotesqr[type],_rotesqr_drv2[type],_rotsize[type]);
}

double SRotDens(double gamma,int type)   // rotational density matrix 
{
   int size = _rotsize[type];

   if (gamma > _rotgrid[type][size-1])   // replace with extrapolation ?  
   return (_rotdens[type][size-1]);

   if (gamma < _rotgrid[type][0])
   { 
//    return (_rotdens[type][0]);
      double rl = _rotgrid [type][0];
      double rr = _rotgrid [type][1];
   
      double salpha = (_rotdens[type][1]    - _rotdens[type][0]   )/(rr-rl);
      double sbeta  = (_rotdens[type][0]*rr - _rotdens[type][1]*rl)/(rr-rl);
     
      return (salpha*gamma + sbeta);
   }

   double rdens;
   splint(_rotgrid[type],_rotdens[type],_rotdens_drv2[type],size,gamma,rdens);
   return rdens;
}

double SRotDensDeriv(double gamma,int type) // the derivs of the rotational density matrix 
{
   int  size = _rotsize[type];

   if (gamma > _rotgrid[type][size-1])      // replace with extrapolation ?  
// return (_rotderv[type][size-1]);
   return 0.0;


   if (gamma < _rotgrid[type][0])
   { 
//    return (_rotderv[type][0]);
      double rl = _rotgrid [type][0];
      double rr = _rotgrid [type][1];
   
      double salpha = (_rotderv[type][1]    - _rotderv[type][0]   )/(rr-rl);
      double sbeta  = (_rotderv[type][0]*rr - _rotderv[type][1]*rl)/(rr-rl);
     
      return (salpha*gamma + sbeta);
   }

   double rderv;
   splint(_rotgrid[type],_rotderv[type],_rotderv_drv2[type],size,gamma,rderv);
   return rderv;
}

double SRotDensEsqrt(double gamma,int type) // the 2nd derivs of the rotational density matrix 
{
   int  size = _rotsize[type];

   if (gamma > _rotgrid[type][size-1])      // replace with extrapolation ?  
// return (_rotesqr[type][size-1]);
   return 0.0;


   if (gamma < _rotgrid[type][0])
   {
//    return (_rotesqr[type][0]);
      double rl = _rotgrid [type][0];
      double rr = _rotgrid [type][1];

      double salpha = (_rotesqr[type][1]    - _rotesqr[type][0]   )/(rr-rl);
      double sbeta  = (_rotesqr[type][0]*rr - _rotesqr[type][1]*rl)/(rr-rl);

      return (salpha*gamma + sbeta);
   }

   double resqr;
   splint(_rotgrid[type],_rotesqr[type],_rotesqr_drv2[type],size,gamma,resqr);
   return resqr;
}

double SPot1D(double r,int atype) 
// it shoud be SPot1D(r,type0,type1) or  SPot1D(r,ind) with ind =type0*NumbTypes+type1
{
   int size = _psize1D[atype];

   if (r >= _pgrid1D[atype][size-1])                //  large distances
   return (-_c6[atype]/pow(r,6.0));
// return 0.0;

   if (r <= _pgrid1D[atype][0])               
   return (_unode[atype]*exp(-_alpha[atype]*r));    //  small distances

   double spot;
   splint(_pgrid1D[atype],_poten1D[atype],_pderiv2[atype],size,r,spot);
   return spot;
}

/*
// revised by Hui Li
double DLPot2D(double r, double cost, int type)
{
   double drmin = _drgrid2D[type][0];
   double dcmin = _dcgrid2D[type][0];

   int drsize   = _drsize2D[type];
   int dcsize   = _dcsize2D[type];

   int  ir = (int)floor((r    - drmin)/_ddelta_r[type]);
   int  ic = (int)floor((cost - dcmin)/_ddelta_c[type]);
  
// CHECK INDEX RANGE  [need to use extrapolation here]
   if (ir<0) ir = 0;
   else if (ir>=(drsize-1)) ir = (drsize - 2); // need to define (ir+1)
	   
   if (ic<0) ic = 0;
   else if (ic>=(dcsize-1)) ic = (dcsize - 2); // need to define (ic+1)
  
// START linear interpolation --------------
  
   double ** dpot = _dpoten2D[type];

   double y1 = dpot[ir][ic];
   double y2 = dpot[ir+1][ic];
   double y3 = dpot[ir+1][ic+1];
   double y4 = dpot[ir][ic+1];

   double r1 = _drgrid2D[type][ir];
   double r2 = _drgrid2D[type][ir+1];

   double c1 = _dcgrid2D[type][ic];
   double c2 = _dcgrid2D[type][ic+1];
   
// CHECK possible division by zero
    
   double dr   = (r   - r1)/(r2-r1);
   double dc   = (cost- c1)/(c2-c1);
   
   double dlpot = (1.0-dr)*(1.0-dc)*y1 + dr*(1.0-dc)*y2 + dr*dc*y3 + (1.0-dr)*dc*y4;

   return dlpot; 
}
*/


double LPot2D(double r, double cost, int type)
{
   double rmin = _rgrid2D[type][0];
   double cmin = _cgrid2D[type][0];

   int rsize   = _rsize2D[type];
   int csize   = _csize2D[type];

   int  ir = (int)floor((r    - rmin)/_delta_r[type]);
   int  ic = (int)floor((cost - cmin)/_delta_c[type]);
  
// CHECK INDEX RANGE  [need to use extrapolation here]
   if (ir<0) ir = 0;
   else if (ir>=(rsize-1)) ir = (rsize - 2); // need to define (ir+1)
	   
   if (ic<0) ic = 0;
   else if (ic>=(csize-1)) ic = (csize - 2); // need to define (ic+1)
  
// START linear interpolation --------------
  
   double ** pot = _poten2D[type];

   double y1 = pot[ir][ic];
   double y2 = pot[ir+1][ic];
   double y3 = pot[ir+1][ic+1];
   double y4 = pot[ir][ic+1];

   double r1 = _rgrid2D[type][ir];
   double r2 = _rgrid2D[type][ir+1];

   double c1 = _cgrid2D[type][ic];
   double c2 = _cgrid2D[type][ic+1];
   
// CHECK possible division by zero
    
   double dr   = (r   - r1)/(r2-r1);
   double dc   = (cost- c1)/(c2-c1);
   
   double lpot = (1.0-dr)*(1.0-dc)*y1 + dr*(1.0-dc)*y2 + dr*dc*y3 + (1.0-dr)*dc*y4;

   return lpot; 
}

int get_filesize(const char fname [])
// count lines with data in input file
{	
   const char *_proc_=__func__;    // "file_size"; 

   ifstream fid(fname,ios::in);   

   if (!fid.good())
   _io_error(_proc_,IO_ERR_FOPEN,fname);

   string stmp;

   int count=0;
   while (fid>>stmp)
   {
      if ((stmp != "") && (stmp != COMMENTS)) count++;
      getline(fid,stmp,'\n');
   } 

   fid.close();
   return (count);

// fid.clear();             // can be used to reset input file 
// fid.seekg(0,ios::beg);
}

void read_datafile(const char fname[],double *grid,double *data)
//  fname:      first   column grid points 
//	        second  column data points
{
   const char *_proc_=__func__;    // "read_DataFile"; 
   
   ifstream fid(fname,ios::in);   

   if (!fid.good())
   _io_error(_proc_,IO_ERR_FOPEN,fname);

   string sgrid;  
   string sdata;
  
   int count=0;
   while (fid>>sgrid)
   {
       if ((sgrid != "") && (sgrid != COMMENTS))  // skip comments and empty lines
       {                                          // should be compatible with get_filesize()
           grid [count] = strtod(sgrid.c_str(),NULL);

           fid>>sdata;
           data [count] = strtod(sdata.c_str(),NULL);

           count++;
       }  
       getline(fid,sdata,'\n');              // skip the rest of the line 
   }
// if (count!= maxsize) 
// nrerror(_proc_,"Wrong size of data file");
   fid.close();
}

void read_datafile(const char fname[],double *grid,double *data0,double *data1)
//  fname:      first   column grid points 
//	        second  column data points
//	        third   column data points
{
   const char *_proc_=__func__;    // "read_DataFile"; 
   
   ifstream fid(fname,ios::in);   

   if (!fid.good())
   _io_error(_proc_,IO_ERR_FOPEN,fname);

   string sgrid;  
   string sdata;
  
   int count=0;
   while (fid>>sgrid)
   {
       if ((sgrid != "") && (sgrid != COMMENTS))  // skip comments and empty lines
       {                                          // should be compatible with get_filesize()
           grid  [count] = strtod(sgrid.c_str(),NULL);

           fid>>sdata;
           data0 [count] = strtod(sdata.c_str(),NULL);

           fid>>sdata;
           data1 [count] = strtod(sdata.c_str(),NULL);

           count++;
       }  
       getline(fid,sdata,'\n');              // skip the rest of the line 
   }
// if (count!= maxsize) 
// nrerror(_proc_,"Wrong size of data file");
   fid.close();
}
void read_datafile(const char fname[],double *grid,double *data0,double *data1,double *data2)
//  fname:      first   column grid points 
//              second  column data points
//              third   column data points
//              fourth  column data points
{
   const char *_proc_=__func__;    // "read_DataFile"; 

   ifstream fid(fname,ios::in);

   if (!fid.good())
   _io_error(_proc_,IO_ERR_FOPEN,fname);

   string sgrid;
   string sdata;

   int count=0;
   while (fid>>sgrid)
   {
       if ((sgrid != "") && (sgrid != COMMENTS))  // skip comments and empty lines
       {                                          // should be compatible with get_filesize()
           grid  [count] = strtod(sgrid.c_str(),NULL);

           fid>>sdata;
           data0 [count] = strtod(sdata.c_str(),NULL);

           fid>>sdata;
           data1 [count] = strtod(sdata.c_str(),NULL);

           fid>>sdata;
           data2 [count] = strtod(sdata.c_str(),NULL);

           count++;
       }
       getline(fid,sdata,'\n');              // skip the rest of the line 
   }
// if (count!= maxsize) 
// nrerror(_proc_,"Wrong size of data file");
   fid.close();
}
