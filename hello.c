//
//         main()
//

#include <stdlib.h>
//#include <stdlib>
#include <sstream>
#include <iomanip>

#include "mc_confg.h"
#include "mc_setup.h"
#include "mc_input.h"
#include "mc_randg.h"
#include "mc_utils.h"
#include "mc_poten.h"
#include "mc_piqmc.h"
#include "mc_qworm.h"
#include "mc_estim.h"

#include <stdio.h>
#include <mpi.h>

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

double _dpot_total;  // potential energy differences, global average  added by Hui Li
double _pot_total;  // kinetic   energy, global average
double _kin_total;  // potential energy, global average 

double _brot;       // rotational kin energy, block average
double _rot_total;  // rotational kin energy, global average

fstream _feng;      // save accumulated energy

//---------------- ESTIMATORS ---------------

void SaveEnergy    (const char [],double,long int); // block average
void SaveSumEnergy (double,double);                 // accumulated average

void InitTotalAverage(void);
void DoneTotalAverage(void);

int main(int argc, char *argvc[])
{

  int numprocs, rank, namelen;
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Get_processor_name(processor_name, &namelen);

  printf("Process %d on %s out of %d\n", rank, processor_name, numprocs);

  MPI_Finalize();
}

