#ifndef __MAIN_H__
#define __MAIN_H__

extern /* readonly */ CProxy_Cell cellArray;
extern /* readonly */ CProxy_Compute computeArray;
extern /* readonly */ CProxy_StaticSchedule staticSch;
extern /* readonly */ int cellArrayDimX;
extern /* readonly */ int cellArrayDimY;
extern /* readonly */ int cellArrayDimZ;

#include <map>
#include <set>
#include <list>
#include <cassert>

struct Obj {
  int objType;
  CkIndex6D idx;

  void pup(PUP::er &p) {
    p | objType;
    p | idx;
  }
};

struct StateNode {
  Obj thisObject, relObject;
  int taskid, taskType, pe, waitTask;

  void pup(PUP::er &p) {
    p | thisObject;
    p | relObject;
    p | taskid;
    p | taskType;
    p | pe;
    p | waitTask;
  }
};

struct SendTo {
  int pe;
  std::list<Obj> sendTo;
};
struct CommInfo {
  int iter, spe, rpe;
  std::list<SendTo> sends;
};

struct OnePerPE : public CkArrayMap {
  OnePerPE() { }
  int procNum(int arrayHdl, const CkArrayIndex &elt) {
    //CkPrintf("procNum = %d\n", *(int*)elt.data());
    assert(*(int*)elt.data() < CkNumPes() &&
           *(int*)elt.data() >= 0);
    return *(int*)elt.data();
  }
};

struct ComputeMap : public CkArrayMap {
  ComputeMap() { }
  int procNum(int arrayHdl, const CkArrayIndex &idx) {
    const short *coor = idx.indexShorts; //(const short*)idx.data();
    int proc = map(coor);
    //CkPrintf("coord = %d %d %d %d %d %d -> proc %d\n",
    //coor[0], coor[1], coor[2], coor[3], coor[4], coor[5], proc);
    return proc;
  }
  int map(const short* coor) {
#include "computeInit"
    else {
      return 0;
    }
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
    else return 99999;
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
