#include "defs.h"
#include "leanmd.decl.h"
#include "Patch.h"
#include "Compute.h"
#include "physics.h"
#include "ckmulticast.h"

extern /* readonly */ CProxy_Main mainProxy;
extern /* readonly */ CProxy_Patch patchArray;
extern /* readonly */ CProxy_Compute computeArray;
extern /* readonly */ CkGroupID mCastGrpID;

extern /* readonly */ int patchArrayDimX;	// Number of Chare Rows
extern /* readonly */ int patchArrayDimY;	// Number of Chare Columns
extern /* readonly */ int patchArrayDimZ;
extern /* readonly */ int finalStepCount; 
extern /* readonly */ BigReal stepTime; 

// Compute - Default constructor
Compute::Compute() {
  cellCount = 0;
  bmsgLenAll = -1;
  stepCount = 0;
  usesAtSync = CmiTrue;
}

Compute::Compute(CkMigrateMessage *msg): CBase_Compute(msg)  { 
  usesAtSync = CmiTrue;
  delete msg;
}

// Entry method to receive vector of particles
void Compute::interact(ParticleDataMsg *msg){
  int i;
  double energy = 0;

  // self interaction check
  if (thisIndex.x1 ==thisIndex.x2 && thisIndex.y1 ==thisIndex.y2 && thisIndex.z1 ==thisIndex.z2) {
    stepCount++;
    bool doatSync = false;
    bmsgLenAll = -1;
    if (msg->doAtSync){
      doatSync = true;
    }
    CkGetSectionInfo(cookie1,msg);
    energy = calcInternalForces(msg, &cookie1, stepCount);
    if(stepCount == 1 || stepCount == finalStepCount)
      contribute(sizeof(double),&energy,CkReduction::sum_double,CkCallback(CkIndex_Main::energySumP(NULL),mainProxy));
    if(doatSync)
      AtSync();
  } else {
    if (cellCount == 0) {
      bufferedMsg = msg;
      bmsgLenAll = bufferedMsg->lengthAll;
      cellCount++;
    } else if (cellCount == 1) {
      // if both particle sets are received, compute interaction
      cellCount = 0;
      stepCount++;
      bool doatSync = false;
      bmsgLenAll = -1;
      if (msg->doAtSync){
        doatSync = true;
      }
      if (bufferedMsg->x*patchArrayDimY*patchArrayDimZ + bufferedMsg->y*patchArrayDimZ + bufferedMsg->z < msg->x*patchArrayDimY*patchArrayDimZ + msg->y*patchArrayDimZ + msg->z){ 
        if (bufferedMsg->lengthAll <= msg->lengthAll)
          energy = calcPairForces(bufferedMsg, msg, &cookie1, &cookie2,stepTime);
        else
          energy = calcPairForces(msg, bufferedMsg, &cookie2, &cookie1,stepTime);
      } else {
        if (bufferedMsg->lengthAll <= msg->lengthAll)
          energy = calcPairForces(bufferedMsg, msg, &cookie2, &cookie1,stepTime);
        else
          energy = calcPairForces(msg, bufferedMsg, &cookie1, &cookie2,stepTime);
      }
      if(stepCount == 1 || stepCount == finalStepCount)
	contribute(sizeof(double),&energy,CkReduction::sum_double,CkCallback(CkIndex_Main::energySumP(NULL),mainProxy));
      bufferedMsg = NULL;
      if(doatSync)
        AtSync();
    }
  }
}

