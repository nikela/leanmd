#ifndef __MAIN_H__
#define __MAIN_H__

extern /* readonly */ CProxy_Cell cellArray;
extern /* readonly */ CProxy_Compute computeArray;

extern /* readonly */ int cellArrayDimX;
extern /* readonly */ int cellArrayDimY;
extern /* readonly */ int cellArrayDimZ;

#include <map>
#include <set>
#include <list>

//central controller chare
class Main : public CBase_Main {
  double startTime;

  private:
    Main_SDAG_CODE

  public:
    Main(CkArgMsg* msg);
    Main(CkMigrateMessage* msg);
    void pup(PUP::er &p);
};
#endif
