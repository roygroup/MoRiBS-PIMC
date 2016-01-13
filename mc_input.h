#ifndef _MC_INPUT_H
#define _MC_INOUT_H 1
#include <stdlib.h>

const char FINPUT[]="qmc.input";   // input params file

// some parameters from this file can be stored in a config file

#include "mc_confg.h"

//enum IOStatus {IORead=0,IOWrite=1}; // not portable ?
const int IORead  = 0;
const int IOWrite = 1; 

void IOReadParams (const char [],int &);  
void StatusIO(int, const char []);       // save/restore status of simulation
void ConfigIO(int, const char []);       // save/restore init configuration
void TablesIO(int, const char []);       // save/restore permutation tables
void QWormsIO(int, const char []);       // save/restore status of the worm

extern string MasterDir;
extern string OutputDir;
extern string FNPrefix;

extern string MCFileName;                // OutputDir+FNPrefix 
 
//----IO errors----------------------------------------------

const char IO_ERR_FOPEN[] = "Can't open input file";
const char IO_ERR_FEXST[] = "File already exists:";
const char IO_ERR_WMODE[] = "Wrong IO mode";
const char IO_ERR_SCALL[] = "Error of the system call";


inline void _io_error(const char proc [],const char error[], const char file_name[])
{
// cerr 
   cout << endl;
   cout << proc << " : " <<error<<"  ["<<file_name<<"]" << endl;
   cout << endl;
   exit(1);
}

// ---- PARAMETERS TO CONTROL INPUT/OUTPUT --------

const int  IO_BLOCKNUMB_WIDTH =  3;
const int  IO_PRECISION       =  6;      
const int  IO_WIDTH           = 14;  // width > precision + 4 + 3    
const int  IO_WIDTH_BLOCK     =  4;  // io width to save counters    

void io_setout(fstream &, int=IO_PRECISION);

int FileExist (const char []);

void IOFileBackUp(const char []);
void IOFileDelete(const char []);

// ------ XYZ DATA FILES --------------------------

const char IO_COM_XYZ[] = "xyz format:  [atom type]  x y z (Angstrom) ";   
const char IO_EXT_XYZ[] = ".xyz";        // extension
       
void IOxyz(int, const char []);     
 
const char IO_EXT_DMP[] = ".dump"; // extra extension for data dumps 

// ----------- ESTIMATORS ------------------------

const char IO_EXT_GRA []    = ".gra";   // radial  pair distribution
const char IO_EXT_GRI []    = ".gri";   // radial  density around impurity
const char IO_EXT_GRT []    = ".grt";   // theta angular density around impurity
const char IO_EXT_GRC []    = ".grc";   // chi angular density around impurity
const char IO_EXT_GTC []    = ".gtc";   // theta and chi angular density around impurity
const char IO_EXT_REP []    = ".eulphi";   // relative euler angle phi
const char IO_EXT_REC []    = ".eulchi";   // relative euler angle chi
const char IO_EXT_RET []    = ".eulthe";   // relative euler angle chi

const char IO_EXT_RCF []    = ".rcf";   // rotational correlation functions

const char IO_EXT_DENS2D [] = ".g2d";   // 2D density
const char IO_EXT_DENS3D [] = ".g3d";   // 3D density

const char IO_EXT_ENG []    = ".eng";   // energy
const char IO_EXT_PRL []    = ".prl";   // distribution of exchange lengths
const char IO_EXT_SUP []    = ".sup";   // area estimators 
const char IO_EXT_MFFSUP3D []    = ".mffs3d";   // area estimators for non-linear rotor in dopant-fixed frame
const char IO_EXT_SFFSUP3D []    = ".sffs3d";   // area estimators for non-linear rotor in space-fixed frame

const char FPERMU []    = "permutation.tab";   // permutation table

const char IO_SUM [] = "_sum";          // file name postfix for accum averages

#endif  //mc_input.h

