#include "time.h"
#include "ckmulticast.h"

#include "defs.h"
#include "leanmd.decl.h"
#include "Main.h"
#include "Patch.h"
#include "Compute.h"

/* readonly */ CProxy_Main mainProxy;
/* readonly */ CProxy_Patch patchArray;
/* readonly */ CProxy_Compute computeArray;
/* readonly */ CkGroupID mCastGrpID;

/* readonly */ int patchArrayDimX;	// Number of Chares in X
/* readonly */ int patchArrayDimY;	// Number of Chares in Y
/* readonly */ int patchArrayDimZ;	// Number of Chares in Z
/* readonly */ int finalStepCount; 
/* readonly */ int firstLdbStep; 
/* readonly */ int ldbPeriod; 
/* readonly */ int ftPeriod; 
/* readonly */ BigReal stepTime; 

// Entry point of Charm++ application
Main::Main(CkArgMsg* msg) {
  LBTurnInstrumentOff();
  stepTime = CmiWallTimer();
  CkPrintf("\nLENNARD JONES MOLECULAR DYNAMICS START UP ...\n");

  patchArrayDimX = PATCHARRAY_DIM_X;
  patchArrayDimY = PATCHARRAY_DIM_Y;
  patchArrayDimZ = PATCHARRAY_DIM_Z;
  finalStepCount = DEFAULT_FINALSTEPCOUNT;
  firstLdbStep = DEFAULT_FIRST_LDB;
  ldbPeriod = DEFAULT_LDB_PERIOD;
  ftPeriod = DEFAULT_FT_PERIOD;

  mainProxy = thisProxy;
  phase = 0;

  int bFactor = 2;

  int numPes = CkNumPes();
  int currPe = -1, pe;

  //TODO Do something about configurable paramters

  // initializing the 3D patch array

  patchArray = CProxy_Patch::ckNew();
  for (int x=0; x<patchArrayDimX; x++)
    for (int y=0; y<patchArrayDimY; y++)
      for (int z=0; z<patchArrayDimZ; z++) {
        pe = (++currPe) % numPes;
        patchArray(x, y, z).insert(pe);
      }
  patchArray.doneInserting();

  CkPrintf("\nNUMBER OF PATCHES: %d X %d X %d .... CREATED\n", patchArrayDimX, patchArrayDimY, patchArrayDimZ);

  mCastGrpID = CProxy_CkMulticastMgr::ckNew(bFactor);

  // initializing the 6D compute array
  computeArray = CProxy_Compute::ckNew();
  for (int x=0; x<patchArrayDimX; x++)
    for (int y=0; y<patchArrayDimY; y++)
      for (int z=0; z<patchArrayDimZ; z++)
	patchArray(x, y, z).createComputes();

  delete msg;
}

// Constructor for chare object migration
Main::Main(CkMigrateMessage* msg): CBase_Main(msg) { }

void Main::lbBarrier(){
  patchArray.resume();
}

void Main::ftBarrier(){
  CkCallback cb(CkIndex_Patch::ftresume(), patchArray);
  CkStartMemCheckpoint(cb);
}

void Main::allDone() {
  CkPrintf("SIMULATION COMPLETE. \n\n");  CkExit();
}

void Main::startUpDone() {
  switch(phase) {
    case 0:
      computeArray.doneInserting();
      CkPrintf("NUMBER OF COMPUTES: %d .... CREATED\n", (NUM_NEIGHBORS/2+1) * patchArrayDimX * patchArrayDimY * patchArrayDimZ);
      phase++;
      patchArray.createSection();
      break;

    case 1:
      CkPrintf("MULTICAST SECTIONS .... CREATED\n");

      CkPrintf("STARTING SIMULATION .... \n\n");
      patchArray.run();
      break;
  }
}

#include "leanmd.def.h"
