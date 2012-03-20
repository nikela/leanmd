#include "defs.h"
#include "leanmd.decl.h"
#include "ComputePME.h"
#include "Cell.h"
#include <vector>

//compute - Default constructor
ComputePME::ComputePME() {
  __sdag_init();
  stepCount = 1;
  energy[0] = energy[1] = 0;
  setMigratable(false);
}

ComputePME::ComputePME(CkMigrateMessage *msg): CBase_ComputePME(msg)  {
  __sdag_init();
  setMigratable(false);
  delete msg;
}

//interaction between two cells
//mcast1 attached to message from lower id cell
void ComputePME::interact(ParticleDataMsg *msg){
  std::vector<vec3> forces(msg->lengthAll);
  cellArray[thisIndex].receivePMEForces(&forces[0], msg->lengthAll);
  delete msg;
}

//pack important information if I am moving
void ComputePME::pup(PUP::er &p) {
  CBase_ComputePME::pup(p);
  __sdag_pup(p);
  p | stepCount;
  PUParray(p, energy, 2);
}
