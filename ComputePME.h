#ifndef LEANMD_COMPUTE_PME_H
#define LEANMD_COMPUTE_PME_H

#include "defs.h"
#include "leanmd.decl.h"

extern /* readonly */ int finalStepCount;
extern /* readonly */ CProxy_Main mainProxy;

//class representing the interaction agents between a couple of cells
class ComputePME : public CBase_ComputePME {
private:
  ComputePME_SDAG_CODE
  int stepCount;  //current step number
  double energy[2]; //store potential energy
  ParticleDataMsg *bufferedMsg; //copy of first message received for interaction

public:
  ComputePME();
  ComputePME(CkMigrateMessage *msg);
  void pup(PUP::er &p);

  void interact(ParticleDataMsg *msg);
};

#endif /* LEANMD_COMPUTE_PME_H */
