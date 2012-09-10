#include "defs.h"
#include "leanmd.decl.h"
#include "Cell.h"
#include "Compute.h"
#include "physics.h"
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

  energyP = calcInternalForces(msg, stepCount);

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

  energyP = calcPairForces(msg, bufferedMsg, stepCount);

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
  PUParray(p, energy, 2);
}
