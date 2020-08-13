#include <string>
#include "time.h"
#include "Main.h"
#include "Cell.h"
#include "Compute.h"

/* readonly */ CProxy_Main mainProxy;
/* readonly */ CProxy_Cell cellArray;
/* readonly */ CProxy_Compute computeArray;

/* readonly */ int cellArrayDimX;
/* readonly */ int cellArrayDimY;
/* readonly */ int cellArrayDimZ;
/* readonly */ int finalStepCount;
/* readonly */ int firstLdbStep;
/* readonly */ int ldbPeriod;
/* readonly */ int checkptFreq;
/* readonly */ int checkptStrategy;
/* readonly */ std::string logs;
/* readonly */ int bFactor;

// Entry point of Charm++ application
Main::Main(CkArgMsg* m) {
  CkPrintf("\nLENNARD JONES MOLECULAR DYNAMICS START UP ...\n");

  // Set variable values to a default set
  cellArrayDimX = CELLARRAY_DIM_X;
  cellArrayDimY = CELLARRAY_DIM_Y;
  cellArrayDimZ = CELLARRAY_DIM_Z;
  finalStepCount = DEFAULT_FINALSTEPCOUNT;
  firstLdbStep = DEFAULT_FIRST_LDB;
  ldbPeriod = DEFAULT_LDB_PERIOD;
  checkptFreq = DEFAULT_FT_PERIOD;
  checkptStrategy = 1;

  mainProxy = thisProxy;

  bFactor = 4;

  int numPes = CkNumPes();
  int currPe = -1, pe;
  int cur_arg = 1;

  CkPrintf("\nInput Parameters ...\n");

  // Read user parameters
  // Number of cells in each dimension
  if (m->argc > cur_arg) {
    cellArrayDimX = atoi(m->argv[cur_arg++]);
    cellArrayDimY = atoi(m->argv[cur_arg++]);
    cellArrayDimZ = atoi(m->argv[cur_arg++]);
    CkPrintf("Cell Array Dimension X: %d, Y: %d, Z: %d of size %d, %d, %d\n",
        cellArrayDimX, cellArrayDimY, cellArrayDimZ, CELL_SIZE_X, CELL_SIZE_Y, CELL_SIZE_Z);
  }

  // Number of steps in simulation
  if (m->argc > cur_arg) {
    finalStepCount = atoi(m->argv[cur_arg++]);
    CkPrintf("Final Step Count: %d\n", finalStepCount);
  }

  // Step after which load balancing starts
  if (m->argc > cur_arg) {
    firstLdbStep = atoi(m->argv[cur_arg++]);
    CkPrintf("First LB Step: %d\n", firstLdbStep);
  }

  // Periodicity of load balancing
  if (m->argc > cur_arg) {
    ldbPeriod = atoi(m->argv[cur_arg++]);
    CkPrintf("LB Period: %d\n", ldbPeriod);
  }

  // Periodicity of checkpointing
  if (m->argc > cur_arg) {
    checkptFreq = atoi(m->argv[cur_arg++]);
    CkPrintf("FT Period: %d\n", checkptFreq);
  }

  // Choose the checkpointing strategy
  if (m->argc > cur_arg) {
  	checkptStrategy = 0;
    logs = m->argv[cur_arg];
  }

  // Create Cell to PE mapping
  CProxy_CellMap cellMap = CProxy_CellMap::ckNew(cellArrayDimX, cellArrayDimY,
      cellArrayDimZ);
  CkArrayOptions opts(cellArrayDimX, cellArrayDimY, cellArrayDimZ);
  opts.setMap(cellMap);

  // Create a 3D Patch array (with a uniform distribution)
  cellArray = CProxy_Cell::ckNew(opts);

  // Create an empty 6D computer array to be filled in by Cells
  computeArray = CProxy_Compute::ckNew();

  cellArray.createComputes();
  CkPrintf("\nCells: %d X %d X %d ... created\n", cellArrayDimX,
    cellArrayDimY, cellArrayDimZ);

  delete m;
}

// Constructor for chare object migration
Main::Main(CkMigrateMessage* msg) : CBase_Main(msg) {}

// PUP routine in case the main chare moves
void Main::pup(PUP::er &p) {
  CBase_Main::pup(p);
  __sdag_pup(p);
}

#include "leanmd.def.h"
