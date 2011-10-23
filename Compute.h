#ifndef __COMPUTE_H__
#define __COMPUTE_H__

#include "defs.h"

//class representing the interaction agents between a couple of cells
class Compute : public CBase_Compute {
  private:
    int cellCount;  // to count the number of interact() calls
    int bmsgLenAll;
    int stepCount;  //current step number
    ParticleDataMsg *bufferedMsg; //copy of first message received for interaction
    //handles to differentiate the two multicast sections I am part of
    CkSectionInfo cookie1;     
    CkSectionInfo cookie2;

  public:
    Compute();
    Compute(CkMigrateMessage *msg);
    void pup(PUP::er &p);

    void interact(ParticleDataMsg *msg);
    void ResumeFromSync() { }           
};

#endif
