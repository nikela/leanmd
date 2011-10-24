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

/* readonly */ int patchArrayDimX;
/* readonly */ int patchArrayDimY;
/* readonly */ int patchArrayDimZ;
/* readonly */ int finalStepCount; 
/* readonly */ int firstLdbStep; 
/* readonly */ int ldbPeriod; 

// Entry point of Charm++ application
Main::Main(CkArgMsg* m) {
   __sdag_init();
  CkPrintf("\nLENNARD JONES MOLECULAR DYNAMICS START UP ...\n");

  //set variable values to a default set
  patchArrayDimX = PATCHARRAY_DIM_X;
  patchArrayDimY = PATCHARRAY_DIM_Y;
  patchArrayDimZ = PATCHARRAY_DIM_Z;
  finalStepCount = DEFAULT_FINALSTEPCOUNT;
  firstLdbStep = DEFAULT_FIRST_LDB;
  ldbPeriod = DEFAULT_LDB_PERIOD;

  mainProxy = thisProxy;
  energy = prevEnergy = 0;
  testFailed = 0;

  //branch factor for spanning tree of multicast
  int bFactor = 4;
  //creating the multicast spanning tree
  mCastGrpID = CProxy_CkMulticastMgr::ckNew(bFactor);

  int numPes = CkNumPes();
  int currPe = -1, pe;
  int cur_arg = 1;

  CkPrintf("\nInput Parameters...\n");

  //read user parameters
  //number of patches/cells in each dimension
  if (m->argc > cur_arg) {
    patchArrayDimX=atoi(m->argv[cur_arg++]);
    patchArrayDimY=atoi(m->argv[cur_arg++]);
    patchArrayDimZ=atoi(m->argv[cur_arg++]);
    CkPrintf("Patch Array Dimension X:%d Y:%d Z:%d\n",patchArrayDimX,patchArrayDimY,patchArrayDimZ);
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

  //initializing the 3D patch array
  patchArray = CProxy_Patch::ckNew();
  for (int x=0; x<patchArrayDimX; x++)
    for (int y=0; y<patchArrayDimY; y++)
      for (int z=0; z<patchArrayDimZ; z++) {
        pe = (++currPe) % numPes;
        patchArray(x, y, z).insert(pe);
      }
  patchArray.doneInserting();

  CkPrintf("\nPatches: %d X %d X %d .... created\n", patchArrayDimX, patchArrayDimY, patchArrayDimZ);

  //initializing the 6D compute array
  computeArray = CProxy_Compute::ckNew();
  for (int x=0; x<patchArrayDimX; x++)
    for (int y=0; y<patchArrayDimY; y++)
      for (int z=0; z<patchArrayDimZ; z++)
        patchArray(x, y, z).createComputes();

  thisProxy.run();
  delete m;
}

//constructor for chare object migration
Main::Main(CkMigrateMessage* msg): CBase_Main(msg) { 
  __sdag_init();  
}

//pup routine incase the main chare moves, pack important information
void Main::pup(PUP::er &p) {
  CBase_Main::pup(p);
   __sdag_pup(p);
  p|energy;
  p|prevEnergy;
  p|testFailed;
}

#include "leanmd.def.h"
