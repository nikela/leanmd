#ifndef __COMPUTE_H__
#define __COMPUTE_H__

#include "defs.h"

//class representing the interaction agents between a couple of cells
class Compute : public CBase_Compute {
  private:
    Compute_SDAG_CODE
    int cellCount;  // to count the number of interact() calls
    int stepCount;  //current step number
    ParticleDataMsg *bufferedMsg; //copy of first message received for interaction
    //handles to differentiate the two multicast sections I am part of
    CkSectionInfo mcast1;     
    CkSectionInfo mcast2;

  public:
    Compute();
    Compute(CkMigrateMessage *msg);
    void pup(PUP::er &p);

    void interact(ParticleDataMsg *msg);
    void ResumeFromSync() { }           
};

#endif
