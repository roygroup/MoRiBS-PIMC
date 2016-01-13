#ifndef _MC_CONFG_H
#define _MC_CONFG_H 1

#include <iostream>
#include <fstream>
#include <string>

using namespace std;

const int MAX_STRING_LENGTH = 80;   // max length of strings in input files 

// --------- PROPERTIES -------------------

const int BOLT  = 0;
const int BOSE  = 1;
const int FERMI = 2;

const char STATISTICS[3][MAX_STRING_LENGTH] = {"BOLTZMANN", "BOSE", "FERMI"};

// --------- SYSTEM -----------------------

const int MAX_NUMBER_TYPES  = 6 ;  // (max+1) number of atom and molecule types
const int MAX_NUMBER_INTER  = 10;  //  max    number of 1D interaction potentials 
const int MAX_NUMBER_ROTDN  = 2;   //  max    number of rotational density matrices  

//--------- FILES  ----------------------

const char COMMENTS[] = "#";        // comment line in output files
const char BLANK[]    = "   ";      // spacing between columns

const char FSTATUS[] = "yw001.stat"; // status of simulations
const char FCONFIG[] = "yw001.conf"; // initial configurations
const char FTABLES[] = "yw001.tabl"; // permutation tables 
const char FRANDOM[] = "yw001.rand"; // status of sprng streams 
const char FQWORMS[] = "yw001.worm"; // status of the worm 

const char FBACKUP_EXTENSION[]=".old";
 
//------- MC-----------------------

const int MC_BLOCK = 0; // block average 
const int MC_TOTAL = 1; // total average (sum over blocks)

//--------  MPI -------------------

const int MPI_MASTER = 0;  // MPI id for the root process 

//----------- ROTATIONS -----------

const int AXIS_X = 0; 
const int AXIS_Y = 1;
const int AXIS_Z = 2; 

const int PHI = 0;      // phi
const int CTH = 1;      // cos(theta)
const int CHI = 2;      // chi

const int NUMB_RCF = 10; // max number Pl(nn) correlation functions

const double RZERO = 1.0e-10;   // zero for the rotational density  matrix
const int    JMAX  = 100;        // used in GetExactRotEnergy(JMAX)
const int    RGRID_SIZE = 1500; // number of grid points for rotational density 

const double PLONE = 0.9999999; // max argument for gsl_sf_legendre_Pl()
//-------------------------------

const int GSLOOP_MAX = 7;  // max length of exchange loops contributing 
                           // to the ground state

const double  HOMEGA = 1.0;

#define DEBUG_WORM  1
#define DEBUG_PIMC  1

#endif  // mc_confg.h
