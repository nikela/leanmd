#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include "mpi-interoperate.h"

int main(int argc, char **argv) {
  int myrank, numranks;
  MPI_Comm newComm;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &numranks);

  CharmLibInit(newComm, argc, argv);

  MPI_Barrier(MPI_COMM_WORLD);
  CharmLibExit();

  printf("called finalize\n");
  MPI_Finalize();
  return 0;
}
