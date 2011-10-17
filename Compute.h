#ifndef __COMPUTE_H__
#define __COMPUTE_H__

#include "defs.h"

// Class representing the interaction agents between a couple of cells
class Compute : public CBase_Compute {
  private:
    int cellCount;  // to count the number of interact() calls
    int bmsgLenAll;
    ParticleDataMsg *bufferedMsg;
    CkSectionInfo cookie1;
    CkSectionInfo cookie2;

  public:
    Compute();
    Compute(CkMigrateMessage *msg);

    void interact(ParticleDataMsg *msg);

    void pup(PUP::er &p) {
      CBase_Compute::pup(p);
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
    void ResumeFromSync();           
};

#endif
