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

class Comm : public CBase_Comm {
  private:
    std::set<int> myCells;
    std::set<int> myComputes;
    std::map< int, std::list< ParticleDataMsg* > > msgs;
    int X, Y, Z, X2, Y2, Z2;
    int linearize6D(CkIndex6D indx) {
      int num = indx.z1*Z2*X*X2*Y*Y2;
      num += indx.z2*X*X2*Y*Y2;
      num += indx.x1*X2*Y*Y2;
      num += indx.x2*Y*Y2;
      num += indx.y1*Y2;
      num += indx.y2;
      return num;
    }

  public:
    Comm() {
      X = CELLARRAY_DIM_X;
      Y = CELLARRAY_DIM_Y;
      Z = CELLARRAY_DIM_Z;
      //X2 =
      //Y2 =
      //Z2 =
    }

    void registerCell(CkIndex3D indx) {
      myCells.insert(indx.z*X*Y+indx.y*X+indx.x);
    }

    void registerCompute(CkIndex6D indx) {
      int num = linearize6D(indx);
      myComputes.insert(num);
    }

    void tryDeliver(ParticleDataMsg* m, CkIndex6D indx) {
      int num = linearize6D(indx);
      if(myComputes.find(num) != myComputes.end())
        deliver(m,indx);
      else
        msgs[num].push_back(m);
    }

    void deliver(ParticleDataMsg* m, CkIndex6D indx) {

    }

    void storeMsg(ParticleDataMsg *m) {

    }

    void checkMsgs(CkIndex6D indx) {
      int num = linearize6D(indx);
      for(std::list< ParticleDataMsg* >::iterator iter = msgs[num].begin();
            iter != msgs[num].end(); ++iter) {
        deliver(*iter, indx);
      }
      msgs[num].clear();
    }

};
#endif
