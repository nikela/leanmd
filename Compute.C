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
extern /* readonly */ double stepTime; 

//compute - Default constructor
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

//entry method to receive vector of particles
void Compute::interact(ParticleDataMsg *msg){
  int i;
  double energy = 0;

  //self interaction check
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
      contribute(sizeof(double),&energy,CkReduction::sum_double,CkCallback(CkReductionTarget(Main,energySumP),mainProxy));
    if(doatSync)
      AtSync();
  } else {
    //check if this is the first message or second
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
      //energy reduction only in begining and end
      if(stepCount == 1 || stepCount == finalStepCount)
        contribute(sizeof(double),&energy,CkReduction::sum_double,CkCallback(CkReductionTarget(Main, energySumP),mainProxy));
      bufferedMsg = NULL;
      if(doatSync)
        AtSync();
    }
  }
}

//pack important information if I am moving
void Compute::pup(PUP::er &p) {
      CBase_Compute::pup(p);
      p | stepCount;
      p | cookie1;
      p | cookie2;
      if (p.isUnpacking() && CkInRestarting()) {
        cookie1.get_redNo() = 0;
        if (!(thisIndex.x1 ==thisIndex.x2 && thisIndex.y1 ==thisIndex.y2 && thisIndex.z1 ==thisIndex.z2))
          cookie2.get_redNo() = 0;
      }
      p | cellCount;
      p | bmsgLenAll;
      int hasMsg = (bmsgLenAll >= 0); // only pup if msg will be used
      p | hasMsg;
      if (hasMsg){
        if (p.isUnpacking())
          bufferedMsg = new (bmsgLenAll) ParticleDataMsg;
        p | *bufferedMsg;
      }
      else
        bufferedMsg = NULL;
}
