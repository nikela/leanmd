#include "defs.h"
#include "leanmd.decl.h"
#include "Cell.h"
#include "Compute.h"
#include "physics.h"
#include "ckmulticast.h"
#include <algorithm>
using std::swap;

extern /* readonly */ CProxy_Main mainProxy;
extern /* readonly */ CkGroupID mCastGrpID;

extern /* readonly */ int cellArrayDimX;
extern /* readonly */ int cellArrayDimY;
extern /* readonly */ int cellArrayDimZ;
extern /* readonly */ int finalStepCount; 

//compute - Default constructor
Compute::Compute() {
  __sdag_init();
  stepCount = 1;
  energy[0] = energy[1] = 0;
  usesAtSync = CmiTrue;
}

Compute::Compute(CkMigrateMessage *msg): CBase_Compute(msg)  { 
  __sdag_init();
  usesAtSync = CmiTrue;
  delete msg;
}

//interaction within a cell
void Compute::selfInteract(ParticleDataMsg *msg){
  double energyP = 0;

  CkGetSectionInfo(mcast1,msg);
  energyP = calcInternalForces(msg, &mcast1, stepCount);

  //energy assignment only in begining and end
  if(stepCount == 1) {
    energy[0] = energyP;
  } else if(stepCount == finalStepCount) {
    energy[1] = energyP;
  }
}

//interaction between two cells
void Compute::interact(ParticleDataMsg *msg){
  double energyP = 0;

  ParticleDataMsg *msgA = msg, *msgB = bufferedMsg;
  CkSectionInfo *handleA = &mcast1, *handleB = &mcast2;
  if (bufferedMsg->x*cellArrayDimY*cellArrayDimZ + bufferedMsg->y*cellArrayDimZ + bufferedMsg->z < msg->x*cellArrayDimY*cellArrayDimZ + msg->y*cellArrayDimZ + msg->z){ 
    swap(handleA, handleB);
  }
  if (bufferedMsg->lengthAll <= msg->lengthAll) {
    swap(msgA, msgB);
    swap(handleA, handleB);
  }
  energyP = calcPairForces(msgA, msgB, handleA, handleB, stepCount);
  bufferedMsg = NULL;

  //energy assignment only in begining and end
  if(stepCount == 1) {
    energy[0] = energyP;
  } else if(stepCount == finalStepCount) {
    energy[1] = energyP;
  }
}

//pack important information if I am moving
void Compute::pup(PUP::er &p) {
  CBase_Compute::pup(p);
  __sdag_pup(p);
  p | stepCount;
  p | mcast1;
  p | mcast2;
  PUParray(p, energy, 2);
  if (p.isUnpacking() && CkInRestarting()) {
    mcast1.get_redNo() = 0;
    if (!(thisIndex.x1 ==thisIndex.x2 && thisIndex.y1 ==thisIndex.y2 && thisIndex.z1 ==thisIndex.z2))
      mcast2.get_redNo() = 0;
  }
  bufferedMsg = NULL;
}
