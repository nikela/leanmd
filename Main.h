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
#include <cassert>

struct OnePerPE : public CkArrayMap {
  OnePerPE() { }
  int procNum(int arrayHdl, const CkArrayIndex &elt) {
    assert(*(int*)elt.data() < CkNumPes() &&
           *(int*)elt.data() >= 0);
    return *(int*)elt.data();
  }
};

struct ComputeMap : public CkArrayMap {
  ComputeMap() { }
  int procNum(int arrayHdl, const CkArrayIndex &idx) {
    const int *coor = idx.data();
    return map(coor);
  }
  int map(const int coor[6]) {
#include "computeInit"
  }
};

struct CellMap : public CkArrayMap {
  CellMap() { }
  int procNum(int arrayHdl, const CkArrayIndex &idx) {
    const int *coor = idx.data();
    return map(coor);
  }
  int map(const int coor[3]) {
#include "cellInit"
  }
};

//central controller chare
class Main : public CBase_Main {
  double startTime;

  private:
    Main_SDAG_CODE

  public:
    Main(CkArgMsg* msg);
    Main(CkMigrateMessage* msg);
    void pup(PUP::er &p);
  void setupFinished();
};
#endif
