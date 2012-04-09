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
  computeTime = CkWallTimer();
  double energyP = 0;

  energyP = calcInternalForces(msg, &mcast1, stepCount);

  //energy assignment only in begining and end
  if(stepCount == 1) {
    energy[0] = energyP;
  } else if(stepCount == finalStepCount) {
    energy[1] = energyP;
  }
  computeTime = CkWallTimer() - computeTime;
}

//interaction between two cells
//mcast1 attached to message from lower id cell
void Compute::interact(ParticleDataMsg *msg){
  computeTime = CkWallTimer();
  double energyP = 0;

  CkSectionInfo *handleA = &mcast1, *handleB = &mcast2;
  if (bufferedMsg->x*cellArrayDimY*cellArrayDimZ + bufferedMsg->y*cellArrayDimZ + bufferedMsg->z < msg->x*cellArrayDimY*cellArrayDimZ + msg->y*cellArrayDimZ + msg->z){ 
    swap(handleA, handleB);
  }
  energyP = calcPairForces(msg, bufferedMsg, handleA, handleB, stepCount);

  //energy assignment only in begining and end
  if(stepCount == 1) {
    energy[0] = energyP;
  } else if(stepCount == finalStepCount) {
    energy[1] = energyP;
  }
  computeTime = CkWallTimer() - computeTime;
}

//pack important information if I am moving
void Compute::pup(PUP::er &p) {
  CBase_Compute::pup(p);
  __sdag_pup(p);
  p | stepCount;
  p | mcast1;
  p | mcast2;
  p | currentState;
  p | stateCount;
  p | thisCompute;
  p | current;
  PUParray(p, energy, 2);
  if (p.isUnpacking() && CkInRestarting()) {
    mcast1.get_redNo() = 0;
    mcast2.get_redNo() = 0;
  }
}
