#include <cmath>
#include <stdio.h>

//Toby set SEED as 985456376 + rank to make different SEEDs for different CPUs
//#define SEED 985456376

#include "sprng.h"	

#include "mc_input.h"
#include "mc_randg.h"
#include "mc_utils.h"
#include "mc_setup.h"
	
const int MAXGENS=15;  // max number of rnd generators [gauss #8,9 N,M #>=10]
int *STREAM[MAXGENS];  // rnd numbers streams


void RandomInit(int mpi_id, int max_mpi_procs)  // sprng init
{            
  int max_streams=MAXGENS*max_mpi_procs;        // number of streams 

  for (int ip=0;ip<MAXGENS;ip++)
  { 
     int streamnum = mpi_id*MAXGENS+ip;
     STREAM[ip] = init_sprng(streamnum,max_streams,SEED,SPRNG_DEFAULT); // initialize stream 
//   print_sprng(stream);	
  }
}

void RandomFree(void) // free memory used by sprng 
{
  for (int ip=0;ip<MAXGENS;ip++)
  free_sprng(STREAM[ip]);     
}

void RandomIO(int tstatus, const char file_name[]) // check point for sprng streams
{
  const char *_proc_=__func__;    // "RandIO()";
  char *fop="r";                  // file operation

  switch (tstatus)
  {
     case IOWrite: fop="w"; break;
     case IORead : fop="r"; break;
     default     :
        nrerror (_proc_,IO_ERR_WMODE); 
        break;
  } 
  
  FILE *fp = fopen(file_name,fop);           
  if(fp == NULL)
  _io_error(_proc_,IO_ERR_FOPEN,file_name);

#ifdef USE_MPI
   not implemented 
#endif

  switch (tstatus)
  {
     case IOWrite: 
         char *bytes;
         for (int ip=0;ip<MAXGENS;ip++)
         {
             int size = pack_sprng(STREAM[ip],&bytes); // pack  a stream state 
             fwrite(&size,1,sizeof(int),fp);           // store # of bytes  
             fwrite(bytes,1,size,fp);                  // store a stream state
         }
         break;
     case IORead:
         char buffer[MAX_PACKED_LENGTH];
         int size;
         for (int ip=0;ip<MAXGENS;ip++)
         {
            fread(&size,1,sizeof(int),fp);
            fread(buffer,1,size,fp);
            STREAM[ip] = unpack_sprng(buffer);
         }
         break;
     default:  
        nrerror (_proc_,IO_ERR_WMODE); 
        break;
  }
  fclose(fp);
}

double rnd(void)
// double random number generator in the range [0,1] 
{
  return sprng(STREAM[0]);
}

double rnd1(void)
// double random number generator in the range [0,1] 
{
  return sprng(STREAM[1]);
}

double rnd2(void)
// double random number generator in the range [0,1] 
{
  return sprng(STREAM[2]);
}

double rnd3(void)
// double random number generator in the range [0,1] 
{
  return sprng(STREAM[3]);
}

double rnd4(void)
// double random number generator in the range [0,1] 
{
  return sprng(STREAM[4]);
}

double rnd5(void)
// double random number generator in the range [0,1] 
{
  return sprng(STREAM[5]);
}

double rnd6(void)
// double random number generator in the range [0,1] 
{
  return sprng(STREAM[6]);
}

double rnd7(void)
// double random number generator in the range [0,1] 
{
  return sprng(STREAM[7]);
}

double gauss(double alpha)
// random numbers according Gaussian distribution exp{-alpha*x^2}
{
    double r1 = sprng(STREAM[8]);
    double r2 = sprng(STREAM[9]);

    double x1 = sqrt(-log(r1))*cos(2.0*M_PI*r2);

//  double x1 = sqrt(-log(r1))*cos(2.0*M_PI*r2);
//  double x2 = sqrt(-log(r1))*sin(2.0*M_PI*r2);

    return (x1/sqrt(alpha));
}

int nrnd1(int n)
// integer random number generator in the range [0,n-1] 
{
  return (int)floor(n*sprng(STREAM[10]));
}

int nrnd2(int n)
// integer random number generator in the range [0,n-1] 
{
  return (int)floor(n*sprng(STREAM[11]));
}

int nrnd3(int n)
// integer random number generator in the range [0,n-1] 
{
  return (int)floor(n*sprng(STREAM[12]));
}

int nrnd4(int n)
// integer random number generator in the range [0,n-1] 
{
  return (int)floor(n*sprng(STREAM[13]));
}
