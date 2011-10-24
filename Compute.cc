#include "defs.h"
#include "leanmd.decl.h"
#include "Patch.h"
#include "Compute.h"
#include "physics.h"
#include "ckmulticast.h"
#include <algorithm>
using std::swap;

extern /* readonly */ CProxy_Main mainProxy;
extern /* readonly */ CkGroupID mCastGrpID;

extern /* readonly */ int patchArrayDimX;
extern /* readonly */ int patchArrayDimY;
extern /* readonly */ int patchArrayDimZ;
extern /* readonly */ int finalStepCount; 

//compute - Default constructor
Compute::Compute() {
  __sdag_init();
  cellCount = 0;
  stepCount = 0;
  usesAtSync = CmiTrue;
}

Compute::Compute(CkMigrateMessage *msg): CBase_Compute(msg)  { 
  __sdag_init();
  usesAtSync = CmiTrue;
  delete msg;
}

//local method to compute forces
void Compute::interact(ParticleDataMsg *msg){
  int i;
  double energy = 0;

  //self interaction check
  if (thisIndex.x1 ==thisIndex.x2 && thisIndex.y1 ==thisIndex.y2 && thisIndex.z1 ==thisIndex.z2) {
    CkGetSectionInfo(mcast1,msg);
    energy = calcInternalForces(msg, &mcast1, stepCount);
    if(stepCount == 0 || stepCount == (finalStepCount-1))
      contribute(sizeof(double),&energy,CkReduction::sum_double,CkCallback(CkReductionTarget(Main,energySumK),mainProxy));
  } else {
    //check if this is the first message or second
    if (cellCount == 0) {
      bufferedMsg = msg;
      cellCount++;
      return;
    }

    // Both particle sets have been received, so compute interaction
    cellCount = 0;

    ParticleDataMsg *msgA = msg, *msgB = bufferedMsg;
    CkSectionInfo *handleA = &mcast1, *handleB = &mcast2;
    if (bufferedMsg->x*patchArrayDimY*patchArrayDimZ + bufferedMsg->y*patchArrayDimZ + bufferedMsg->z < msg->x*patchArrayDimY*patchArrayDimZ + msg->y*patchArrayDimZ + msg->z){ 
      swap(handleA, handleB);
    }
    if (bufferedMsg->lengthAll <= msg->lengthAll) {
      swap(msgA, msgB);
      swap(handleA, handleB);
    }

    energy = calcPairForces(msgA, msgB, handleA, handleB, stepCount);

    //energy reduction only in begining and end
    if(stepCount == 0 || stepCount == (finalStepCount-1))
      contribute(sizeof(double),&energy,CkReduction::sum_double,CkCallback(CkReductionTarget(Main, energySumK),mainProxy));
    bufferedMsg = NULL;
  }
}

//pack important information if I am moving
void Compute::pup(PUP::er &p) {
  CBase_Compute::pup(p);
  __sdag_pup(p);
  p | stepCount;
  p | mcast1;
  p | mcast2;
  if (p.isUnpacking() && CkInRestarting()) {
    mcast1.get_redNo() = 0;
    if (!(thisIndex.x1 ==thisIndex.x2 && thisIndex.y1 ==thisIndex.y2 && thisIndex.z1 ==thisIndex.z2))
      mcast2.get_redNo() = 0;
  }
  p | cellCount;
  bufferedMsg = NULL;
}
