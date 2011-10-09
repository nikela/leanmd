#include "time.h"
#include "physics.h"
#include "mol3d.decl.h"
#include "Main.h"
#include "Patch.h"
#include "Compute.h"

/* readonly */ CProxy_Main mainProxy;
/* readonly */ CProxy_Patch patchArray;
/* readonly */ CProxy_Compute computeArray;
/* readonly */ CkGroupID mCastGrpID;

/* readonly */ int numParts;
/* readonly */ int patchArrayDimX;	// Number of Chares in X
/* readonly */ int patchArrayDimY;	// Number of Chares in Y
/* readonly */ int patchArrayDimZ;	// Number of Chares in Z
/* readonly */ int patchSizeX;
/* readonly */ int patchSizeY;
/* readonly */ int patchSizeZ;
/* readonly */ int ptpCutOff;
/* readonly */ int patchMargin;
/* readonly */ int patchOriginX;
/* readonly */ int patchOriginY;
/* readonly */ int patchOriginZ;
/* readonly */ int migrateStepCount;
/* readonly */ int finalStepCount; 
/* readonly */ int firstLdbStep; 
/* readonly */ int ldbPeriod; 
/* readonly */ int ftPeriod; 
/* readonly */ BigReal stepTime; 
/* readonly */ BigReal timeDelta;
/* readonly */ bool usePairLists;
/* readonly */ bool twoAwayX;
/* readonly */ bool twoAwayY;
/* readonly */ bool twoAwayZ;
/* readonly */ int numNbrs;
/* readonly */ int nbrsX;
/* readonly */ int nbrsY;
/* readonly */ int nbrsZ;

/* readonly */ double A = 1.60694452*pow(10.0, -134);// Force Calculation parameter 1
/* readonly */ double B = 1.031093844*pow(10.0, -77);// Force Calculation parameter 2

// Entry point of Charm++ application
Main::Main(CkArgMsg* msg) {
  LBTurnInstrumentOff();
  stepTime = CmiWallTimer();
  CkPrintf("\nLENNARD JONES MOLECULAR DYNAMICS START UP ...\n");
  numParts = 0;

  usePairLists = USE_PAIRLISTS;
  patchArrayDimX = PATCHARRAY_DIM_X;
  patchArrayDimY = PATCHARRAY_DIM_Y;
  patchArrayDimZ = PATCHARRAY_DIM_Z;
  ptpCutOff = PTP_CUT_OFF;
  patchMargin = PATCH_MARGIN;
  patchSizeX = PATCH_SIZE_X;
  patchSizeY = PATCH_SIZE_Y;
  patchSizeZ = PATCH_SIZE_Z;
  patchOriginX = PATCH_ORIGIN_X;
  patchOriginY = PATCH_ORIGIN_Y;
  patchOriginZ = PATCH_ORIGIN_Z;
  migrateStepCount = MIGRATE_STEPCOUNT;
  finalStepCount = DEFAULT_FINALSTEPCOUNT;
  firstLdbStep = DEFAULT_FIRST_LDB;
  ldbPeriod = DEFAULT_LDB_PERIOD;
  ftPeriod = DEFAULT_FT_PERIOD;
  timeDelta = DEFAULT_DELTA;
  numNbrs = NUM_NEIGHBORS;
  nbrsX = NBRS_X;
  nbrsY = NBRS_Y;
  nbrsZ = NBRS_Z;
  twoAwayX = false;
  twoAwayY = false;
  twoAwayZ = false;

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

#ifdef USE_SECTION_MULTICAST
  mCastGrpID = CProxy_CkMulticastMgr::ckNew(bFactor);
#endif

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
      CkPrintf("NUMBER OF COMPUTES: %d .... CREATED\n", (numNbrs/2+1) * patchArrayDimX * patchArrayDimY * patchArrayDimZ);
#ifdef USE_SECTION_MULTICAST
      phase++;
      patchArray.createSection();
      break;

    case 1:
      CkPrintf("MULTICAST SECTIONS .... CREATED\n");
#endif

      CkPrintf("STARTING SIMULATION .... \n\n");
      patchArray.start();
      break;
  }
}

#include "mol3d.def.h"
