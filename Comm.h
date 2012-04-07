#ifndef __COMM_H__
#define __COMM_H__

#include <cassert>

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
      X2 = X;
      Y2 = Y;
      Z2 = Z;
    }

    void registerCell(CkIndex3D indx) {
      myCells.insert(indx.z*X*Y+indx.y*X+indx.x);
      //checkMsgs(indx);
    }

    void registerCompute(CkIndex6D indx) {
      int num = linearize6D(indx);
      myComputes.insert(num);
      checkMsgs(indx);
    }

    void unregister(CkIndex3D indx) {
      assert(myCells.find(indx.z*X*Y+indx.y*X+indx.x) != myCells.end());
      myCells.erase(indx.z*X*Y+indx.y*X+indx.x);
    }

    void unregister(CkIndex6D indx) {
      int num = linearize6D(indx);
      assert(myComputes.find(num) != myComputes.end());
      myComputes.erase(num);
    }

    void tryDeliver(ParticleDataMsg* m) {
      CkIndex3D cellIndx;
      cellIndx.x = m->x;
      cellIndx.y = m->y;
      cellIndx.z = m->z;
      CkIndex6D indx;//= secGrp[cellIndx][stepCount].ckLocal()->getID();
      int num = linearize6D(indx);
      if(myComputes.find(num) != myComputes.end())
        deliver(m,indx);
      else
        msgs[num].push_back(m);
    }

    void deliver(ParticleDataMsg* m, CkIndex6D indx) {
      computeArray[indx].ckLocal()->interact(m);
      //TODO: When to do a "selfinteract"?
    }

    void sendMsg(ParticleDataMsg *m, CkIndex3D cellIndx) {
      //secGrp[cellIndx][stepCount].ckLocal()->getSection().tryDeliver(m);
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

#endif /*__COMM_H__*/
