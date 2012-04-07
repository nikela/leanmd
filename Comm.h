#ifndef __COMM_H__
#define __COMM_H__

#include <cassert>

class Comm : public CBase_Comm {
  private:
    std::set<int> myCells;
    std::set<int> myComputes;
    std::map< int, std::list< ParticleDataMsg* > > msgs;
    std::map< int, std::list< std::pair< vec3*, int> > > vecmsgs;

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
      checkMsgs(indx);
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

    void tryDeliver(vec3* forces, int n) {
      int num = forces[n-1].x;
      CkIndex3D indx;
      indx.z = num % (X*Y);
      indx.y = (num - indx.z*X*Y) % X;
      indx.x = num - indx.z - indx.y;
      if(myCells.find(num) != myCells.end()) {
        deliver(forces, n, indx);
      }
      else {
        vecmsgs[num].push_back(*(new std::pair<vec3*,int> (forces,n)));
      }
    }

    void tryDeliver(ParticleDataMsg* m) {
      CkIndex3D cellIndx;
      cellIndx.x = m->x;
      cellIndx.y = m->y;
      cellIndx.z = m->z;
      std::list<CkIndex6D> indxs;//= secGrp[cellIndx][stepCount].ckLocal()->getID();

      for (std::list<CkIndex6D>::iterator iter = indxs.begin(); iter != indxs.end(); ++iter) {
        int num = linearize6D(*iter);
        if(myComputes.find(num) != myComputes.end()) {
          deliver(m,*iter);
        }
        else {
          msgs[num].push_back(m);
        }
      }
    }

    void deliver(vec3* forces, int n, CkIndex3D indx) {
      //Remember that this sends n-1 because the nth element is
      //storing the Cell index!
      cellArray[indx].ckLocal()->reduceForces(forces, n-1);
    }

    void deliver(ParticleDataMsg* m, CkIndex6D indx) {
      computeArray[indx].ckLocal()->calculateForces(m);
    }

    void reduceForces(vec3* forces, int n) {
      tryDeliver(forces, n);
    }

    void sendParticles(ParticleDataMsg *m, CkIndex3D cellIndx) {
      //secGrp[cellIndx][stepCount].ckLocal()->getSection().tryDeliver(m);
    }

    void checkMsgs(CkIndex3D indx) {
      int num = indx.z*X*Y+indx.y*X+indx.x;
      for(std::list< std::pair<vec3*,int> >::iterator iter =
            vecmsgs[num].begin(); iter != vecmsgs[num].end(); ++iter) {
        deliver((*iter).first, (*iter).second, indx);
      }
      vecmsgs[num].clear();
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
