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
/* readonly */ double stepTime; 

// Entry point of Charm++ application
Main::Main(CkArgMsg* m) {
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
  energy = prevEnergy = 0;
  testFailed = 0;

  int bFactor = 4;

  int numPes = CkNumPes();
  int currPe = -1, pe;
  int cur_arg = 1;

  CkPrintf("\nInput Parameters...\n");

  if (m->argc > cur_arg) {
    patchArrayDimX=atoi(m->argv[cur_arg++]);
    patchArrayDimY=atoi(m->argv[cur_arg++]);
    patchArrayDimZ=atoi(m->argv[cur_arg++]);
    CkPrintf("Patch Array Dimension X:%d Y:%d Z:%d\n",patchArrayDimX,patchArrayDimY,patchArrayDimZ);
  }

  if (m->argc > cur_arg) {
    finalStepCount=atoi(m->argv[cur_arg++]);
    CkPrintf("Final Step Count:%d\n",finalStepCount);
  }

  if (m->argc > cur_arg) {
    firstLdbStep=atoi(m->argv[cur_arg++]);
    CkPrintf("First LB Step:%d\n",firstLdbStep);
  }

  if (m->argc > cur_arg) {
    ldbPeriod=atoi(m->argv[cur_arg++]);
    CkPrintf("LB Period:%d\n",ldbPeriod);
  }

  if (m->argc > cur_arg) {
    ftPeriod=atoi(m->argv[cur_arg++]);
    CkPrintf("FT Period:%d\n",ldbPeriod);
  }

  // initializing the 3D patch array
  patchArray = CProxy_Patch::ckNew();
  for (int x=0; x<patchArrayDimX; x++)
    for (int y=0; y<patchArrayDimY; y++)
      for (int z=0; z<patchArrayDimZ; z++) {
        pe = (++currPe) % numPes;
        patchArray(x, y, z).insert(pe);
      }
  patchArray.doneInserting();

  CkPrintf("\nPatches: %d X %d X %d .... created\n", patchArrayDimX, patchArrayDimY, patchArrayDimZ);

  mCastGrpID = CProxy_CkMulticastMgr::ckNew(bFactor);

  // initializing the 6D compute array
  computeArray = CProxy_Compute::ckNew();
  for (int x=0; x<patchArrayDimX; x++)
    for (int y=0; y<patchArrayDimY; y++)
      for (int z=0; z<patchArrayDimZ; z++)
	patchArray(x, y, z).createComputes();

  delete m;
}

// Constructor for chare object migration
Main::Main(CkMigrateMessage* msg): CBase_Main(msg) { }

void Main::ftBarrier(){
  CkCallback cb(CkIndex_Patch::ftresume(), patchArray);
  CkStartMemCheckpoint(cb);
}

void Main::allDone() {
  if(testFailed)
    CkPrintf("\nEnergy conservation test failed for maximum allowed variation of %E units.\nSIMULATION UNSUCCESSFULL\n",ENERGY_VAR);  
  else
    CkPrintf("\nEnergy conservation test passed for maximum allowed variation of %E units.\nSIMULATION SUCCESSFULL \n",ENERGY_VAR);
  CkExit();
}

void Main::startUpDone() {
  switch(phase) {
    case 0:
      computeArray.doneInserting();
      CkPrintf("Computes: %d .... created\n", (NUM_NEIGHBORS/2+1) * patchArrayDimX * patchArrayDimY * patchArrayDimZ);
      phase++;
      patchArray.createSection();
      break;

    case 1:
      CkPrintf("Multicast sections .... created\n");

      CkPrintf("Starting simulation .... \n\n");
      patchArray.nextStep();
      break;
  }
}

void Main::energySumK(double energyK) {
  if(energy == 0) {
    energy = energyK;
  } else {
    energy += energyK;
    if(prevEnergy == 0) 
      prevEnergy = energy;
    if(abs(energy-prevEnergy)>ENERGY_VAR) {
      CkPrintf("Energy value has varied significantly from %E to %E\n",prevEnergy,energy);
      testFailed = 1;
    }
    prevEnergy = energy;
    energy = 0;
    //patchArray.testDone(1);
  }
}

void Main::energySumP(double energyP) {
  if(energy == 0) {
    energy = energyP;
  } else {
    energy += energyP;
    if(prevEnergy == 0) 
      prevEnergy = energy;
    if(abs(energy-prevEnergy)>ENERGY_VAR) {
      CkPrintf("Energy value has varied significantly from %E to %E\n",prevEnergy,energy);
      testFailed = 1;
    }
    prevEnergy = energy;
    energy = 0;
    //patchArray.testDone(1);
  }
}

#include "leanmd.def.h"
