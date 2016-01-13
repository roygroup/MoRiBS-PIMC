#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <sstream>
#include <iomanip>

#include "mc_confg.h"
#include "mc_setup.h"
#include "mc_input.h"
#include "mc_utils.h"
#include "mc_poten.h"
#include "mc_qworm.h"

//----- INPUT PARAMETERS----------------------------

// FILE NAMES

const char IO_MASTERDIR[]      = "MASTERDIR";               
const char IO_OUTPUTDIR[]      = "OUTPUTDIR";                 
const char IO_FILENAMEPREFIX[] = "FILENAMEPREFIX";     

// SYSTEM

const char IO_DIMENSION[]      = "DIMENSION";      

const char IO_RESTART[]        = "RESTART";
const char IO_TEMPERATURE[]    = "TEMPERATURE";
const char IO_ATOM[]           = "ATOM";
const char IO_MOLECULE[]       = "MOLECULE";
// Toby's value
const char IO_NONLINEAR[]      = "NONLINEAR";
const char IO_DENSITY[]        = "DENSITY";
const char IO_READMCCOORDS[]   = "READMCCOORDS";

// MC 
const char IO_NUMBEROFSLICES[] = "NUMBEROFSLICES";
const char IO_NUMBEROFPASSES[] = "NUMBEROFPASSES";
const char IO_NUMBEROFBLOCKS[] = "NUMBEROFBLOCKS";
const char IO_INITIALBLOCK[]   = "INITIALBLOCK";

const char IO_REFLECTY[]       = "REFLECTY";
const char IO_REFLECTX[]       = "REFLECTX";
const char IO_REFLECTZ[]       = "REFLECTZ";
const char IO_ROTSYM[]         = "ROTSYM";

const char IO_ROTATION[]       = "ROTATION";
const char IO_ROTDENSI[]       = "ROTDENSI";

const char IO_WORM[]           = "WORM";
const char IO_MINIMAGE[]       = "MINIMAGE";

const char IO_MCSKIP_RATIO[]   = "MCSKIP_RATIO";
const char IO_MCSKIP_TOTAL[]   = "MCSKIP_TOTAL"; 
const char IO_MCSKIP_AVERG[]   = "MCSKIP_AVERG";  

//------ MC STATUS----------------------------------

const char STATUS_STARTBLOCK[] = "STARTBLOCK";

//--------------------------------------------------

string MasterDir;
string OutputDir;
string FNPrefix;

string MCFileName;     // mc output file name (no extension)
 
void IOReadParams(const char in_file[],int & mc_status)
//   mc_status =  0 - new run, 1 - restart  
{
   const char *_proc_=__func__;    // "MCReadInput()"; 

//-- Open file with input parameters ----------------

   ifstream inf(in_file,ios::in);
   if (!inf.good())
   _io_error(_proc_,IO_ERR_FOPEN,in_file);

//-- Read input parameters --------------------------

   int numb   = 0;         // total number of particles in the system
   int type   = 0;         // total number of atom types

   IMPURITY = false;
   WORM     = false;
   MINIMAGE = false;

   ROTATION = false;

   NUMB_ATOMS = 0;         // number of atoms     (structureless particles)
   NUMB_MOLCS = 0;         // number of molecules (rotational degrees of freedom)

   NUMB_ATOMTYPES = 0;     // number of atoms types
   NUMB_MOLCTYPES = 0;     // number of molecules types

   mc_status  = 0;         // 0 - new run, 1 - restart  

   int _wtypes = 0;        // number of atom's types to apply the worm algorithm
   int _rtypes = 0;        // number of molecule's types to sample rot degrees of freedom

   double  _rot_step;      // step for sampling the rotational degrees of freedom
   string _srot_type;      // the molecule type to sample the rotational degrees of freedom
   string _srot_dens;      // the file name with the rotational density matrix

   InitMCCoords = 0;

   MCStartBlock = 0;
   
   string params;
   while (inf>>params)
   {  
     if (params==IO_RESTART)
     {
       mc_status = 1;  // restart 
     }
     else
     if (params==IO_MASTERDIR)
     {
        inf>>MasterDir;
     }
     else
     if (params==IO_OUTPUTDIR)
     {
        inf>>OutputDir;
     }
     else
     if (params==IO_FILENAMEPREFIX)
     {
        inf>>FNPrefix;
     }
     else
     if (params==IO_TEMPERATURE)            
     { 
        inf>>Temperature;
     }
     else 
     if (params==IO_DENSITY)            
     { 
        inf>>Density;
     }
     else 
     if ((params==IO_ATOM)||(params==IO_MOLECULE)||(params==IO_NONLINEAR))
     {
        inf>>MCAtom[type].type;          // [1] 
        inf>>MCAtom[type].numb;          // [2]

        if(MCAtom[type].numb < 0)
        {
           ISPHER = 1;
           MCAtom[type].numb = -MCAtom[type].numb;
        }

        string sstat;
        int     stat; 
        inf>>sstat;                      // [3]
 
        if      (sstat == STATISTICS[BOSE ]) stat = BOSE; 
        else if (sstat == STATISTICS[BOLT ]) stat = BOLT; 
//      else if (sstat == STATISTICS[FERMI]) stat = FERMI; 
        else 
        nrerror(_proc_,"Unknown statistics");

        MCAtom[type].stat = stat;

        inf>>MCAtom[type].mcstep;        // [4]
        inf>>MCAtom[type].levels;        // [5]
        inf>>MCAtom[type].fpot;          // [6]

        string smod;                     // model of interaction 
        int    pmod; 
        inf>>smod;                       // [7]

        if      (smod == PMODEL[PRIMITIVE]) pmod = PRIMITIVE; 
        else if (smod == PMODEL[EFFECTIVE]) pmod = EFFECTIVE; 
        else 
        nrerror(_proc_,"Unknown model of interaction");

        MCAtom[type].pmod = pmod;

// ------- set atom/molecule flag ------------------------- 

        MCAtom[type].molecule = 0;  // atom

        if (params == IO_MOLECULE)
        MCAtom[type].molecule = 1;  // molecular impurity
        else if (params == IO_NONLINEAR)
        MCAtom[type].molecule = 2;
        else                        // atom
        if (IMPURITY) nrerror(_proc_,"Molecules should follow atoms in input file");
//      the latter is important, for example, for the density estimators

        if (MCAtom[type].numb > 0)   // ignore this atom/molecule type if N <= 0 
        {
           if ((MCAtom[type].molecule == 1)||(MCAtom[type].molecule == 2))
           {  
              IMPURITY         = true;
              NUMB_MOLCS      += MCAtom[type].numb;
              NUMB_MOLCTYPES  ++;
           }
           else
           {
              NUMB_ATOMS      += MCAtom[type].numb;
              NUMB_ATOMTYPES  ++;
           }

           numb += MCAtom[type].numb;        
           type++;                         // SHOULD BE THE LAST LINE IN THIS SECTION
        } 
     }
     else 
     if (params==IO_NUMBEROFSLICES)
     {
         inf>>NumbTimes;
     }  
     else 
     if (params==IO_NUMBEROFPASSES)
     {
         inf>>NumberOfMCPasses;
     }  
     else 
     if (params==IO_NUMBEROFBLOCKS)
     {
         inf>>NumberOfMCBlocks>>NumberOfEQBlocks;
     }
     else
     if (params==IO_INITIALBLOCK)
     {
         inf>>MCStartBlock;
     }
     else
     if (params==IO_REFLECTY)
     {
         inf>>IREFLY;
     }
     else
     if (params==IO_REFLECTX)
     {
         inf>>IREFLX;
     }
     if (params==IO_REFLECTZ)
     {
         inf>>IREFLZ;
     }
     else
     if (params==IO_ROTSYM)
     {
         IROTSYM=1;
         inf>>NFOLD_ROT;
     }
     else
     if (params==IO_DIMENSION)
     {
         inf>>NDIM;
     }
     else
     if (params==IO_ROTATION)
     {
        inf>>_srot_type;   //   [1] the molecule type for rotational sampling
 
        inf>>_rot_step;    //   [2] the rotational MC step(s) 
//      inf>>_srot_dens;   //   [3] the file name with the rotational density

        inf>> NumbRotTimes;//   [3]  number of rotational time slices     
 
        ROTATION = true;

       _rtypes ++;
     }
     else
     if (params==IO_ROTDENSI)
     {
        inf>>RotDenType;
        inf>>RotOdEvn;
        inf>>RotEoff;
        inf>>X_Rot;
        inf>>Y_Rot;
        inf>>Z_Rot;
        inf>>RNratio;
     }
     else
     if (params==IO_WORM)
     {
        inf>>Worm.stype;      //     [1]
 
        inf>>Worm.c;          //     [2] 
        inf>>Worm.m;          //     [3]

        WORM = true;

       _wtypes ++;
     }
     else 
     if (params==IO_MINIMAGE)
     {
        MINIMAGE = true;
     }
     else
     if (params==IO_READMCCOORDS)
     {
        InitMCCoords = 1;
     }
     else 
     if (params==IO_MCSKIP_RATIO)
     {
        inf >> MCSKIP_RATIO;
     } 
     else 
     if (params==IO_MCSKIP_TOTAL)
     {
        inf >> MCSKIP_TOTAL;
     } 
     else 
     if (params==IO_MCSKIP_AVERG)
     {
        inf >> MCSKIP_AVERG;
     } 
     else
     {}

     getline(inf,params,'\n');  // skip comments at the end of the line 
   }

   inf.close(); 

// POST READ INITIALIZATION

   NumbAtoms = numb;  // the total number of atoms/molecules
   NumbTypes = type;  // the total number of atoms/molecules types

   int gatom  = 0;    // offset = NumbTimes*gatom
   for (int type=0;type<NumbTypes;type++)
   {
      MCAtom[type].gatom  = gatom;             // the other parts of the code relies on this  
      MCAtom[type].offset = gatom*NumbTimes;   // MCCoords[][offset+NumbTimes*atom+it]            
 
      gatom += MCAtom[type].numb; 
   }

   MaxnTimes                = gatom*NumbTimes;  // (NumbTimes*NumbAtoms)
   MCAtom[NumbTypes].offset = MaxnTimes;        //  used to define MCatom[type+1]

   MCFileName = (OutputDir+FNPrefix);

   if (ROTATION)
   {  
       if (_rtypes != 1)   
       nrerror(_proc_,"Rotational degrees of freedom for one molecule's type only");

       bool found   = false;

       for (int type=0;type<NumbTypes;type++)
       if  (_srot_type == MCAtom[type].type)
       {
          if (MCAtom[type].molecule == 0)
          nrerror(_proc_,"Rotational degrees of freedom for molecules only");

          MCAtom[type].rtstep =  _rot_step;
//        strcpy(MCAtom[type].rdens,_srot_dens.c_str());

          found  = true;  
          break;
       } 
 
       if (!found)
       nrerror(_proc_,"Can't find a particle type to sample rotational degrees of freedom");
   }

   if (Worm.m >=NumbTimes)
   {
      nrerror(_proc_,"Worm algorithm: m should be smaller then M");
   }

   if ((NUMB_ATOMTYPES >1 ) || (NUMB_MOLCTYPES > 1 ))
   nrerror(_proc_,"No more then one atom/molecule type: densities and potential energy");

   if (WORM && (_wtypes != 1))   
   nrerror(_proc_,"Worm algorithm only for one atom type");

#ifdef HOSC_TEST
   if (NumbTypes !=1)
   nrerror(_proc_,"HOSC_TEST: only one atom type");
#endif 

// begin  DUMP -----------------------------------------

   cout << endl;
   cout << "MasterDir   "       <<MasterDir<<endl;
   cout << "OutPut   "          <<OutputDir<<endl;
   cout << "File Name Prefix   "<<FNPrefix<<endl;

   cout << "Temperature "<<Temperature<<endl;
   cout << "Density  "   <<Density<<endl;
   cout << "DIM = "        <<NDIM<<endl;

   int w = 6;
   
   for (int type=0;type<NumbTypes;type++)
   {
      cout << setw(w) << MCAtom[type].type << BLANK;
      cout << setw(w) << MCAtom[type].numb << BLANK;

      cout << setw(10) << STATISTICS[MCAtom[type].stat] << BLANK;

      cout << setw(w) << MCAtom[type].mcstep << BLANK;
      cout << setw(w) << MCAtom[type].levels << BLANK;
      cout << setw(w) << MCAtom[type].fpot   << BLANK;

      cout << setw(w) << PMODEL[MCAtom[type].pmod] << BLANK;

      cout << endl; 
   }

   cout << "Total number of atoms  " << NumbAtoms << endl;
   cout << "Total number of types  " << NumbTypes << endl;

   if (WORM)
   {  
      cout << "WORM"  << BLANK  << Worm.stype;
      cout << setw(w) << Worm.c << BLANK;
      cout << setw(w) << Worm.m << BLANK;
      cout << endl;
   }

   if (ROTATION)
   { 
      cout << endl;
      cout << "ROTATION" <<endl;

      cout << "Rotational step:           "    <<  _rot_step   << endl;
//    cout << "Rotational density matrix: "    << _srot_dens   << endl;
      cout << "Number of Rotational Slices = " << NumbRotTimes << endl;
      cout << "Rotational Density Type = " << RotDenType << endl;
      cout << "RotOdEvn = "<< RotOdEvn <<endl;
      cout << "Rotational Energy Estimator offset = "<< RotEoff <<"CM-1"<<endl;
      cout << "Rotational Constants (CM-1) :"<<X_Rot<<" "<<Y_Rot<<" "<<Z_Rot<<endl;
      cout << "RS / Noya: "<<RNratio<<endl;

      cout << endl;
   }

   if(ROTATION && ISPHER == 1)
   nrerror(_proc_,"ISPHER = 1 is not compatible with ROTATION");

   cout << "NONLINEAR DOPANT SPHERICAL TREATMENT: "<<ISPHER<<endl;

   cout << "Number of Slices = " << NumbTimes << endl;
   cout << "Number of Passes = " << NumberOfMCPasses << endl;
   cout << "Number of Blocks = " << NumberOfMCBlocks << BLANK << NumberOfEQBlocks << endl;
   cout << "Initial Block No.= " << MCStartBlock<<endl;

   cout << endl;
   cout << endl;

   cout << "Number of steps to skip to save ACCEPT RATIO" << BLANK << MCSKIP_RATIO << endl;
   cout << "Number of steps to skip to save ACCUML AVERG" << BLANK << MCSKIP_TOTAL << endl;
   cout << "Number of steps to skip to evaluate AVERAGES" << BLANK << MCSKIP_AVERG << endl;

   cout << endl;
   cout << endl;

   if(IREFLY == 1)
   cout <<"REFLECT WITH RESPECT TO XZ PLANE OF DOPANT WHEN EVALUATING PROPERTIES"<<endl;

   if(IREFLX == 1)
   cout <<"REFLECT WITH RESPECT TO YZ PLANE OF DOPANT WHEN EVALUATING PROPERTIES"<<endl;

   if(IREFLZ == 1)
   cout <<"REFLECT WITH RESPECT TO XY PLANE OF DOPANT WHEN EVALUATING PROPERTIES"<<endl;

   if(IROTSYM ==1)
   cout <<"ROTATIONAL SYMMETRY OF THE DOPANT WITH NFOLD="<<NFOLD_ROT<<endl;

#ifdef HOSC_TEST
   cout << "   <<<<<<  Harmonic oscillator test [HOSC_TEST 1] >>>>>>>   " << endl;
#endif

#ifdef ROTS_TEST
   cout << "   <<<<<<  Free rotor test [ROTS_TEST 1]  >>>>>>>   " << endl;
#endif

   cout << endl;
   cout << endl;

// end    DUMP -----------------------------------------

// return mc_status;
}

void StatusIO(int tstatus, const char file_name[]) 
// read/wrire status of MC run
{
   const char *_proc_=__func__;    // "MCStatus()"; 

   ios::openmode mode;

   switch (tstatus)
   {
      case IOWrite: mode=ios::out;   break;
      case IORead : mode=ios::in;    break;
      default     :
      nrerror (_proc_,IO_ERR_WMODE); break;      
   } 

   fstream fid(file_name,mode);

   if (!fid.good())
   _io_error(_proc_,IO_ERR_FOPEN,file_name);

   string status;
   switch (tstatus)
   {
      case IOWrite: 
         fid<<STATUS_STARTBLOCK<<" "<<MCStartBlock<<endl;
         break;
      case IORead: 
         while (fid>>status)
         {  
            if (status==STATUS_STARTBLOCK)
            fid>>MCStartBlock;
 
            getline(fid,status,'\n');  // skip comments 
         }
         break;
      default :  
         nrerror (_proc_,IO_ERR_WMODE);
         break;      
   } 
  
   fid.close();
}

void ConfigIO(int tstatus, const char file_name[]) 
// read/wrire initial configuration 
{
   const char *_proc_=__func__;    // "MCStatus()"; 

   ios::openmode mode;

   switch (tstatus)
   {
      case IOWrite: mode = ios::out; break;
      case IORead : mode = ios::in;  break;
      default     :
      nrerror (_proc_,IO_ERR_WMODE); break;      
   } 

   fstream fid(file_name,mode);
//   fstream fid(file_name,mode | ios::binary);

   if (!fid.good())
   _io_error(_proc_,IO_ERR_FOPEN,file_name);

   io_setout(fid);  //added by Hui Li

   streamsize size;         
   switch (tstatus)
   {
      case IOWrite: 
         size=sizeof(double)*NumbAtoms*NumbTimes;
         fid.write((char *)&size,sizeof(streamsize));
         fid.write((char *)MCCoords[0],size);
         fid.write((char *)MCCosine[0],size);
//       Toby replaces the above line by
//       fid.write((char *)MCAngles[0],size);
//       to store the three Euler angles
         break;
      case IORead:
         fid.read((char *)&size,sizeof(streamsize));
         fid.read((char *)MCCoords[0],size);
         fid.read((char *)MCCosine[0],size);
//       Toby replaces the above line by
//       fid.read((char *)MCAngles[0],size);
//       to read the three Euler angles
         break;
      default :  
         nrerror (_proc_,IO_ERR_WMODE);
         break;      
   } 
  
   fid.close();
}

void TablesIO(int tstatus, const char file_name[]) 
// read/wrire permutation tables 
{
   const char *_proc_=__func__;    // "TablesIO"; 

   ios::openmode mode;

   switch (tstatus)
   {
      case IOWrite: mode = ios::out; break;
      case IORead : mode = ios::in;  break;
      default     :
      nrerror (_proc_,IO_ERR_WMODE); break;      
   } 

   fstream fid(file_name,mode | ios::binary);

   if (!fid.good())
   _io_error(_proc_,IO_ERR_FOPEN,file_name);

   streamsize size;         
   switch (tstatus)
   {
      case IOWrite: 
         size=sizeof(int)*NumbAtoms;
         fid.write((char *)&size,sizeof(streamsize));
         fid.write((char *)PIndex,size);
         fid.write((char *)RIndex,size);
         break;
      case IORead:
         fid.read((char *)&size,sizeof(streamsize));
         fid.read((char *)PIndex,size);
         fid.read((char *)RIndex,size);
         break;
      default :  
         nrerror (_proc_,IO_ERR_WMODE);
         break;      
   } 
  
   fid.close();
}

void QWormsIO(int tstatus, const char file_name[]) 
// read/write a status of the worm 
{
   const char *_proc_=__func__;  // "QWormsIO()"; 

   ios::openmode mode;

   switch (tstatus)
   {
      case IOWrite: mode = ios::out; break;
      case IORead : mode = ios::in;  break;
      default     :
      nrerror (_proc_,IO_ERR_WMODE); break;      
   } 

   fstream fid(file_name,mode | ios::binary);

   if (!fid.good())
   _io_error(_proc_,IO_ERR_FOPEN,file_name);

// streamsize size;         
   switch (tstatus)
   {
      case IOWrite: 
         fid.write((char *)&Worm,sizeof(TPathWorm));
         break;
      case IORead:
         fid.read((char *)&Worm,sizeof(TPathWorm));
         break;
      default :  
         nrerror (_proc_,IO_ERR_WMODE);
         break;      
   } 
  
   fid.close();
}

void IOxyz(int tstatus, const char file_name[]) 
{
   const char *_proc_=__func__;    // "IOxyz"; 

//---------------- Open  ------------

   ios::openmode mode;

   switch (tstatus)
   {
      case IOWrite: mode = ios::out;  break;
      case IORead : mode = ios::in;   break;
      default     :
      nrerror (_proc_,IO_ERR_WMODE);  break;      
   } 

   fstream fid(file_name,mode);

   if (!fid.good())
   _io_error(_proc_,IO_ERR_FOPEN,file_name);

   io_setout(fid);

//---------------- Read/Write ------------

   stringstream stype; 
   string       sbuff;    
  
   int offset;

   int type = 0;
   int atom = 0;  // first atom # will be 1, NOT 0 
   switch (tstatus)
   {
      case IOWrite:
         fid<<MaxnTimes<<endl;                    // total number of "atoms" 
         fid<<COMMENTS<<BLANK<<IO_COM_XYZ<<endl;  // comments

         for (int it=0;it<MaxnTimes;it++)
         {
            if (it==MCAtom[type+1].offset) {type++; atom=0;}    // new atom type
            if ((it-MCAtom[type].offset)%NumbTimes==0) atom++;  // new atom

            stype.str(""); stype<<MCAtom[type].type<<atom;      // atom label
            fid<<setw(5)<<stype.str()<<BLANK; 
                                                 
            for (int id=0;id<NDIM;id++)
            {
            fid<<setw(IO_WIDTH)<<MCCoords[id][it]<<BLANK;
            fid<<setw(IO_WIDTH)<<MCCosine[id][it]<<BLANK;
//          Toby replaces the above line by
//          fid<<setw(IO_WIDTH)<<MCAngles[id][it]<<BLANK;
//          by doing that, Toby stores the three Euler angles, not the unit vector of the two angles orientation
            }
            fid<<endl;
         } 

/*  the same as above, but explicit sum over atoms and types

         for (int type=0;type<NumbTypes;type++)
         for (int atom=0;atom<MCAtom[type].numb;atom++)
         for (int it=0;it<NumbTimes;it++)
         {
            stype.str(""); stype<<MCAtom[type].type<<atom;     // atom label
            fid<<setw(5)<<stype.str()<<BLANK; 
           
            offset=MCAtom[type].offset+NumbTimes*atom;
 
            for (int id=0;id<NDIM;id++) 
            fid<<MCCoords[id][offset+it]<<BLANK;
            fid<<endl;
         } 
*/ 
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

void IOFileBackUp(const char file_name[])
{
#ifdef DEBUG_PIMC
   const char *_proc_=__func__;    
#endif 

   string cmd="cp";
 
   if (FileExist(file_name))
   {    
      cmd += " " + (string)file_name;
      cmd += " " + (string)file_name+(string)FBACKUP_EXTENSION;
      int status = system(cmd.c_str());

#ifdef DEBUG_PIMC
      if (status < 0)
      nrerror(_proc_,IO_ERR_SCALL);
#endif
   } 
}

void IOFileDelete(const char file_name[])
{          
#ifdef DEBUG_PIMC
   const char *_proc_=__func__;    
#endif
 
   string cmd="rm";
 
   if (FileExist(file_name))
   {    
      cmd += " " + (string)file_name;
      int status = system(cmd.c_str());

#ifdef DEBUG_PIMC
      if (status < 0)
      nrerror(_proc_,IO_ERR_SCALL);
#endif
   } 
}

int FileExist (const char fileName [])
{
   FILE *infile = fopen (fileName, "r");
// ifstream  infile(fileName,ios::in); 
 
   int ret_code=0;
   if (infile == NULL)   
       ret_code = 0;    // file doesn't exist
   else 
   {
       ret_code = 1;     // file exists
       fclose (infile);
   }
   return (ret_code);
}

void io_setout(fstream &out, int set_precision)
{ 
   out<<setprecision(set_precision);
   out<<setiosflags (ios::scientific);
// out<<setiosflags (ios::right);
}
