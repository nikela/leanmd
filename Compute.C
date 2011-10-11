#include "defs.h"
#include "leanmd.decl.h"
#include "Patch.h"
#include "Compute.h"
#include "physics.h"
#ifdef USE_SECTION_MULTICAST
#include "ckmulticast.h"
#endif

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
  LBTurnInstrumentOff();
  cellCount = 0;
  numLists = -1;
  bmsgLenAll = -1;
  usesAtSync = CmiTrue;
}

Compute::Compute(CkMigrateMessage *msg): CBase_Compute(msg)  { 
  usesAtSync = CmiTrue;
  delete msg;
}
  
// Entry method to receive vector of particles
void Compute::interact(ParticleDataMsg *msg){
  int i;

  // self interaction check
  if (thisIndex.x1 ==thisIndex.x2 && thisIndex.y1 ==thisIndex.y2 && thisIndex.z1 ==thisIndex.z2) {
    bool doatSync = false;
    bmsgLenAll = -1;
    if (msg->doAtSync){
      doatSync = true;
    }
    if (msg->lbOn)
      LBTurnInstrumentOn();
    CkGetSectionInfo(cookie1,msg);
    calcInternalForces(msg, &cookie1);
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
      bool doatSync = false;
      bmsgLenAll = -1;
      if (msg->doAtSync){
	doatSync = true;
      }
      if (msg->lbOn)
	LBTurnInstrumentOn();
      if (USE_PAIRLISTS){
	  if (bufferedMsg->x*patchArrayDimY*patchArrayDimZ + bufferedMsg->y*patchArrayDimZ + bufferedMsg->z < msg->x*patchArrayDimY*patchArrayDimZ + msg->y*patchArrayDimZ + msg->z){ 
	    if (bufferedMsg->lengthAll <= msg->lengthAll)
	      pairList = calcPairForcesPL(bufferedMsg, msg, pairList, &numLists, &cookie1, &cookie2);
	    else
	      pairList = calcPairForcesPL(msg, bufferedMsg, pairList, &numLists, &cookie2, &cookie1);
	  }
	  else{
	    if (bufferedMsg->lengthAll <= msg->lengthAll)
	      pairList = calcPairForcesPL(bufferedMsg, msg, pairList, &numLists, &cookie2, &cookie1);
	    else
	      pairList = calcPairForcesPL(msg, bufferedMsg, pairList, &numLists, &cookie1, &cookie2);
	  }
      }
      else {
	if (bufferedMsg->x*patchArrayDimY*patchArrayDimZ + bufferedMsg->y*patchArrayDimZ + bufferedMsg->z < msg->x*patchArrayDimY*patchArrayDimZ + msg->y*patchArrayDimZ + msg->z){ 
	  if (bufferedMsg->lengthAll <= msg->lengthAll)
	    calcPairForces(bufferedMsg, msg, &cookie1, &cookie2);
	  else
	    calcPairForces(msg, bufferedMsg, &cookie2, &cookie1);
	}
	else{
	  if (bufferedMsg->lengthAll <= msg->lengthAll)
	    calcPairForces(bufferedMsg, msg, &cookie2, &cookie1);
	  else
	    calcPairForces(msg, bufferedMsg, &cookie1, &cookie2);
	}
      }
      bufferedMsg = NULL;
      if(doatSync)
	AtSync();
    }
  }
}


void Compute::ResumeFromSync(){
  LBTurnInstrumentOff();
}
