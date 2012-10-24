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
  stepCount = 1;
  energy[0] = energy[1] = 0;
  usesAtSync = CmiTrue;
}

Compute::Compute(CkMigrateMessage *msg): CBase_Compute(msg)  { 
  usesAtSync = CmiTrue;
  delete msg;
}

//interaction within a cell
void Compute::selfInteract(ParticleDataMsg *msg){
  double energyP = 0;
  vec3 * fmsg;

  energyP = calcInternalForces(msg,  stepCount, fmsg);

  //energy assignment only in begining and end
  if(stepCount == 1) {
    energy[0] = energyP;
  } else if(stepCount == finalStepCount) {
    energy[1] = energyP;
  }

  //contribute to force reduction
  CkMulticastMgr *mCastGrp = CProxy_CkMulticastMgr(mCastGrpID).ckLocalBranch();
  CkGetSectionInfo(mcast1, msg);
  mCastGrp->contribute(sizeof(vec3)*msg->lengthAll, fmsg, CkReduction::sum_double, mcast1);

  delete msg;
  delete [] fmsg;
}

//interaction between two cells
void Compute::interact(ParticleDataMsg *msg1, ParticleDataMsg *msg2){
  double energyP = 0;
  vec3 *fmsg1, *fmsg2;

  CkSectionInfo *handleA = &mcast1, *handleB = &mcast2;
  if (msg2->x*cellArrayDimY*cellArrayDimZ + msg2->y*cellArrayDimZ + msg2->z < msg1->x*cellArrayDimY*cellArrayDimZ + msg1->y*cellArrayDimZ + msg1->z){ 
    swap(handleA, handleB);
  }
  energyP = calcPairForces(msg1, msg2, stepCount, fmsg1, fmsg2);

  //energy assignment only in begining and end
  if(stepCount == 1) {
    energy[0] = energyP;
  } else if(stepCount == finalStepCount) {
    energy[1] = energyP;
  }
  
  //contribute to force reduction
  CkMulticastMgr *mCastGrp = CProxy_CkMulticastMgr(mCastGrpID).ckLocalBranch();
  CkGetSectionInfo(*handleA, msg1);
  mCastGrp->contribute(sizeof(vec3)*msg1->lengthAll, fmsg1, CkReduction::sum_double, *handleA);
  CkGetSectionInfo(*handleB, msg2);
  mCastGrp->contribute(sizeof(vec3)*msg2->lengthAll, fmsg2, CkReduction::sum_double, *handleB);

  delete msg1;
  delete msg2;
  delete [] fmsg1;
  delete [] fmsg2;
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
    mcast2.get_redNo() = 0;
  }
}
