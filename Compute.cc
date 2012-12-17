#include "defs.h"
#include "Cell.h"
#include "Compute.h"
#include "physics.h"
#include "ckmulticast.h"
#include <algorithm>
using std::swap;

//compute - Default constructor
Compute::Compute() : stepCount(1) {
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
  std::vector<vec3> force1;

  energyP = calcInternalForces(msg, stepCount, force1);

  //energy assignment only in begining and end
  if (stepCount == 1) energy[0] = energyP;
  else if (stepCount == finalStepCount) energy[1] = energyP;

  //contribute to force reduction
  CkMulticastMgr *mCastGrp = CProxy_CkMulticastMgr(mCastGrpID).ckLocalBranch();
  CkGetSectionInfo(mcast1, msg);
  mCastGrp->contribute(sizeof(vec3)*msg->lengthAll, &force1[0], CkReduction::sum_double, mcast1);

  delete msg;
}

//interaction between two cells
void Compute::interact(ParticleDataMsg *msg1, ParticleDataMsg *msg2){
  CkSectionInfo *handleA = &mcast1, *handleB = &mcast2;
  if (msg2->x * cellArrayDimY * cellArrayDimZ + msg2->y * cellArrayDimZ + msg2->z <
      msg1->x * cellArrayDimY * cellArrayDimZ + msg1->y * cellArrayDimZ + msg1->z)
    swap(handleA, handleB);

  std::vector<vec3> force1, force2;
  double energyP = calcPairForces(msg1, msg2, stepCount, force1, force2);

  //energy assignment only in begining and end
  if (stepCount == 1) energy[0] = energyP;
  else if (stepCount == finalStepCount) energy[1] = energyP;

  //contribute to force reduction
  CkMulticastMgr *mCastGrp = CProxy_CkMulticastMgr(mCastGrpID).ckLocalBranch();
  CkGetSectionInfo(*handleA, msg1);
  mCastGrp->contribute(sizeof(vec3)*msg1->lengthAll, &force1[0], CkReduction::sum_double, *handleA);
  CkGetSectionInfo(*handleB, msg2);
  mCastGrp->contribute(sizeof(vec3)*msg2->lengthAll, &force2[0], CkReduction::sum_double, *handleB);

  delete msg1;
  delete msg2;
}

void Compute::registerResumeClient(){
  	CkCallback _cb(CkIndex_Compute::resumeFromChkp(),thisProxy(thisIndex.x1,thisIndex.y1,thisIndex.z1,thisIndex.x2,thisIndex.y2,thisIndex.z2));
	setChkpResumeClient(_cb);
}

//pack important information if I am moving
void Compute::pup(PUP::er &p) {
  
//  if(p.isChecking())
//  	CkPrintf("[%d][%d] pup compute\n",CmiMyPartition(),CkMyPe());

//  if(p.isChecking())
//	  p.skip();
  CBase_Compute::pup(p);
  __sdag_pup(p);
  
//  if(p.isChecking())
//	  p.resume();

  
  p | stepCount;
  
//  if(p.isChecking())
//  	CkPrintf("[%d][%d] pup stepCount %d\n",CmiMyPartition(),CkMyPe(),stepCount);
  
  PUParray(p, energy, 2);
  
//  if(p.isChecking())
//  	CkPrintf("[%d][%d] pup compute ends\n",CmiMyPartition(),CkMyPe());
  
//  if(p.isChecking())
//	  p.skip();
  p | mcast1;
  p | mcast2;
  
  if (p.isUnpacking() && CkInRestarting()) {
    mcast1.set_redNo(0); mcast2.set_redNo(0);
	registerResumeClient();
  }
  
  
 // if(p.isChecking())
//	  p.resume();
}
