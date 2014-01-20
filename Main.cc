#include <string>
#include "time.h"
#include "ckmulticast.h"
#include "Main.h"
#include "Cell.h"
#include "Compute.h"

/* readonly */ CProxy_Main mainProxy;
/* readonly */ CProxy_Cell cellArray;
/* readonly */ CProxy_Compute computeArray;
/* readonly */ CkGroupID mCastGrpID;

/* readonly */ int cellArrayDimX;
/* readonly */ int cellArrayDimY;
/* readonly */ int cellArrayDimZ;
/* readonly */ int finalStepCount; 
/* readonly */ int firstLdbStep; 
/* readonly */ int ldbPeriod;
/* readonly */ int checkptFreq; 
/* readonly */ int checkptStrategy;
/* readonly */ std::string logs;

/*  readonly*/ int KAWAY_X;
/*  readonly*/ int KAWAY_Y;
/*  readonly*/ int KAWAY_Z;
/*  readonly*/ int NBRS_X;
/*  readonly*/ int NBRS_Y;
/*  readonly*/ int NBRS_Z;
/*  readonly*/ int NUM_NEIGHBORS;
/*  readonly*/ int  CELL_SIZE_X;
/*  readonly*/ int  CELL_SIZE_Y;
/*  readonly*/ int  CELL_SIZE_Z;



// Entry point of Charm++ application
Main::Main(CkArgMsg* m) {
  int problemDimX, problemDimY, problemDimZ;
  CkPrintf("\nLENNARD JONES MOLECULAR DYNAMICS START UP ...\n");

  //set variable values to a default set
  problemDimX = cellArrayDimX = CELLARRAY_DIM_X;
  problemDimY = cellArrayDimY = CELLARRAY_DIM_Y;
  problemDimZ = cellArrayDimZ = CELLARRAY_DIM_Z;
  finalStepCount = DEFAULT_FINALSTEPCOUNT;
  firstLdbStep = DEFAULT_FIRST_LDB;
  ldbPeriod = DEFAULT_LDB_PERIOD;
  checkptFreq = DEFAULT_FT_PERIOD;

  mainProxy = thisProxy;

  //branch factor for spanning tree of multicast
  int bFactor = 4;

  int numPes = CkNumPes();
  int currPe = -1, pe;
  int cur_arg = 1;

  CkPrintf("\nInput Parameters...\n");

  //read user parameters
  //number of cells in each dimension
  if (m->argc > cur_arg) {
    problemDimX = atoi(m->argv[cur_arg++]);
    problemDimY = atoi(m->argv[cur_arg++]);
    problemDimZ = atoi(m->argv[cur_arg++]);
  }

  //number of steps in simulation
  if (m->argc > cur_arg) {
    finalStepCount=atoi(m->argv[cur_arg++]);
    CkPrintf("Final Step Count:%d\n",finalStepCount);
  }

  //step after which load balancing starts
  if (m->argc > cur_arg) {
    firstLdbStep=atoi(m->argv[cur_arg++]);
    CkPrintf("First LB Step:%d\n",firstLdbStep);
  }

  //periodicity of load balancing
  if (m->argc > cur_arg) {
    ldbPeriod=atoi(m->argv[cur_arg++]);
    CkPrintf("LB Period:%d\n",ldbPeriod);
  }

  //periodicity of checkpointing
  if (m->argc > cur_arg) {
    checkptFreq=atoi(m->argv[cur_arg++]);
    CkPrintf("FT Period:%d\n",checkptFreq);
  }

  checkptStrategy = 1;
  //choose the checkpointing strategy use in disk checkpointing
  if (m->argc > cur_arg) {
  	checkptStrategy = 0;
    logs = m->argv[cur_arg++];
  }
    KAWAY_X = 1;
    KAWAY_Y = 1;
    KAWAY_Z = 1;
  //decomposition , Away X
  if (m->argc > cur_arg) {
    KAWAY_X = atoi(m->argv[cur_arg++]);
    CkPrintf("Away X : %d ", KAWAY_X);
  }

  if (m->argc > cur_arg) {
    KAWAY_Y = atoi(m->argv[cur_arg++]);
    CkPrintf("Away Y : %d ", KAWAY_Y);
  }

  if (m->argc > cur_arg) {
    KAWAY_Z = atoi(m->argv[cur_arg++]);
    CkPrintf("Away Z : %d ", KAWAY_Z);
  }
  
  if (m->argc > cur_arg) {
    bFactor = atoi(m->argv[cur_arg]);
    CkPrintf("Branch factor: %d ", bFactor);
  }
  //creating the multicast spanning tree
  mCastGrpID = CProxy_CkMulticastMgr::ckNew(bFactor);

  NBRS_X = (2*KAWAY_X+1);
  NBRS_Y = (2*KAWAY_Y+1);
  NBRS_Z = (2*KAWAY_Z+1);
  NUM_NEIGHBORS = (NBRS_X * NBRS_Y * NBRS_Z);
  CELL_SIZE_X = (PTP_CUT_OFF + CELL_MARGIN)/KAWAY_X;
  CELL_SIZE_Y = (PTP_CUT_OFF + CELL_MARGIN)/KAWAY_Y;
  CELL_SIZE_Z = (PTP_CUT_OFF + CELL_MARGIN)/KAWAY_Z;

  cellArrayDimX = problemDimX * KAWAY_X;
  cellArrayDimY = problemDimY * KAWAY_Y;
  cellArrayDimZ = problemDimZ * KAWAY_Z;
  CkPrintf("Cell Array Dimension X:%d Y:%d Z:%d of size %d %d %d\n",cellArrayDimX,cellArrayDimY,cellArrayDimZ,CELL_SIZE_X,CELL_SIZE_Y,CELL_SIZE_Z);
  delete m;
  cellArray = CProxy_Cell::ckNew();
  computeArray = CProxy_Compute::ckNew();
  thisProxy.init();
}

void Main::init()
{
  //initializing the 3D Patch array (with a uniform distribution) and 6D compute array
  int patchCount = 0;
  float ratio = ((float)CkNumPes() - 1)/(cellArrayDimX*cellArrayDimY*cellArrayDimZ);
  for (int x=0; x<cellArrayDimX; x++)
    for (int y=0; y<cellArrayDimY; y++)
      for (int z=0; z<cellArrayDimZ; z++) {
        cellArray(x, y, z).insert((int)(patchCount++ * ratio));
        cellArray(x, y, z).createComputes();
      }

  cellArray.doneInserting();
  CkPrintf("\nCells: %d X %d X %d .... created\n", cellArrayDimX, cellArrayDimY, cellArrayDimZ);

}

//constructor for chare object migration
Main::Main(CkMigrateMessage* msg): CBase_Main(msg) { 
}

//pup routine incase the main chare moves, pack important information
void Main::pup(PUP::er &p) {
  CBase_Main::pup(p);
  __sdag_pup(p);
}

#include "leanmd.def.h"
