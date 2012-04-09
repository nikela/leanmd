#ifndef __COMM_H__
#define __COMM_H__

extern /*readonly*/ CProxy_Comm commProxy;
extern /* readonly */ int cellArrayDimX;
extern /* readonly */ int cellArrayDimY;
extern /* readonly */ int cellArrayDimZ;

#include "defs.h"
#include "leanmd.decl.h"
#include "Compute.h"
#include <map>
#include <set>
#include <list>
#include <cassert>

extern /* readonly */ CProxy_Compute computeArray;

class Comm : public CBase_Comm {
  private:
    std::set<int> myCells;
    std::set<int> myComputes;
    std::map< int, std::list< ParticleDataMsg* > > msgs;
    std::map< int, std::list< std::pair< vec3*, int> > > vecmsgs;

  public:
    int X, Y, Z, X2, Y2, Z2;
    int linearize6D(CkIndex6D indx) {
      int
      num  = indx.z1 * Z2 * X  * X2 * Y  * Y2;
      num += indx.z2 * X  * X2 * Y  * Y2;
      num += indx.x1 * X2 * Y  * Y2;
      num += indx.x2 * Y  * Y2;
      num += indx.y1 * Y2;
      num += indx.y2;
      return num;
    }

    Comm() {
      X = CELLARRAY_DIM_X;
      Y = CELLARRAY_DIM_Y;
      Z = CELLARRAY_DIM_Z;
      X2 = X;
      Y2 = Y;
      Z2 = Z;
    }

    Comm(CkMigrateMessage*) { }

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

    void calculateForces(ParticleDataMsg *m) {
      CkIndex3D cellIndx;
      cellIndx.x = m->x;
      cellIndx.y = m->y;
      cellIndx.z = m->z;
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

  void warmup() {

  }
};

#include <cassert>

//# object types [ 0 -> CELL,
//#                1 -> COMPUTE,
//#                2 -> NONE ]

//# task types   [ 0 -> CELL_SEND,
//#                1 -> CELL_UPDATE,
//#                2 -> COMPUTE,
//#                3 -> MIGRATE ]

static int linearize3D(CkIndex3D idx) {
  return idx.x * cellArrayDimY * cellArrayDimZ +
         idx.y * cellArrayDimZ +
         idx.z;
}

struct Obj {
  int objType;
  CkIndex6D idx;
};

struct StateNode {
  Obj thisObject, relObject;
  int taskid, taskType, pe, waitTask;
};

struct SendTo {
  int pe;
  std::list<Obj> sendTo;
};
struct CommInfo {
  int iter, spe, rpe;
  std::list<SendTo> sends;
};

struct StaticSchedule : public CBase_StaticSchedule {
  std::map<int, std::list<StateNode> > cellTransition;
  std::map<int, std::list<StateNode> > computeTransition;

  //# iter <cell-obj> (spe,rpe) -> [ (pe,[obj,iter]) ]
  std::map<int, std::map<int, CommInfo> > commMap;

  StaticSchedule() {
    if (CkMyPe() == 0) {
      CkPrintf("reading schedule\n");
      fflush(stdout);
    }

    int curVal = 0, x;

    // read from file;
    FILE* schedule = fopen("schedule.dag", "r");
    StateNode node;
    while (fscanf(schedule, "%d ", &x) == 1) {
      curVal++;
      if (x == -1) {
        //printf("resetting curVal = %d, objtype = %d\n", curVal, node.thisObject.objType);
	curVal = 0;
      } else {
        
	switch (curVal) {
	case 1: node.taskid = x; break;
	case 2: node.thisObject.objType = x; break;
        case 3: node.thisObject.idx.x1 = x; break;
        case 4: node.thisObject.idx.x2 = x; break;
        case 5: node.thisObject.idx.y1 = x; break;
        case 6: node.thisObject.idx.y2 = x; break;
        case 7: node.thisObject.idx.z1 = x; break;
        case 8: node.thisObject.idx.z2 = x; break;
        case 9: node.taskType = x; break;
        case 10: node.pe = x; break;
        case 11: node.waitTask = x; break;
	case 12: node.relObject.objType = x; break;
        case 13: node.relObject.idx.x1 = x; break;
        case 14: node.relObject.idx.x2 = x; break;
        case 15: node.relObject.idx.y1 = x; break;
        case 16: node.relObject.idx.y2 = x; break;
        case 17: node.relObject.idx.z1 = x; break;
        case 18: node.relObject.idx.z2 = x;
          //printf("curVal = %d, objtype = %d\n", curVal, node.thisObject.objType);
          if (node.thisObject.objType == 0) {
            //printf("found cell\n");
            // cell
            cellTransition[commProxy[CkMyPe()].ckLocal()->linearize6D(node.thisObject.idx)].push_back(node);
          } else if (node.thisObject.objType == 1) {
            //printf("found compute\n");
            // compute
            computeTransition[commProxy[CkMyPe()].ckLocal()->linearize6D(node.thisObject.idx)].push_back(node);
          } else {
            assert(0);
          }
	  curVal = 0;
	  break;
	}
      }
    }
    fclose(schedule);

    curVal = 0;

    CommInfo cinfo;
    CkIndex3D sender;
    int peSendCount = 0;

    printf("reading collMap\n");

    FILE* procMap = fopen("collMap", "r");

    while (fscanf(procMap, "%d ", &x) == 1) {
      curVal++;
      if (x == -1) {
	curVal = 0;
      } else {
	switch (curVal) {
	case 1: cinfo.iter = x; break;
	case 2: assert(x == 0);  break; // is the CELL_TYPE
	case 3: sender.x = x; break;
        case 4: break;
        case 5: sender.y = x; break;
        case 6: break;
        case 7: sender.z = x; break;
        case 8: break;
        case 9: cinfo.spe = x; break;
        case 10: cinfo.rpe = x; break;
        case 11: peSendCount = x;
          SendTo send;
          int curVal2 = 0;
          while (fscanf(procMap, "%d ", &x) == 1) {
            curVal2++;
            int objCount = 0;
            switch (curVal2) {
            case 1: send.pe = x; break;
            case 2: objCount = x;

              Obj obj;
              int curVal3 = 0;
              while (fscanf(procMap, "%d ", &x) == 1) {
                curVal3++;
                switch (curVal3) {
                case 1: obj.objType = x; break;
                case 2: obj.idx.x1 = x; break;
                case 3: obj.idx.x2 = x; break;
                case 4: obj.idx.y1 = x; break;
                case 5: obj.idx.y2 = x; break;
                case 6: obj.idx.z1 = x; break;
                case 7: obj.idx.z2 = x;
                  send.sendTo.push_back(obj);
                  curVal3 = 0;
                  objCount--;
                  break;
                }
                if (objCount == 0) break;
              }

              //printf("inserting into cinfo a send, size of sendTo is %d\n",
              //send.sendTo.size());

              cinfo.sends.push_back(send);
              curVal2 = 0;
              peSendCount--;
              break;
            }

            if (peSendCount == 0) break;
          }

          //printf("inserting into commMap, iter = %d, sender = %d\n",
          //cinfo.iter, linearize3D(sender));
          commMap[cinfo.iter][linearize3D(sender)] = cinfo;
          cinfo.sends.clear();
          curVal = 0;

          break;
        }
      }
    }

    fclose(procMap);

    printf("finished reading collMap\n");

    if (CkMyPe() == 0) {
      CkPrintf("finished reading schedule\n");
      fflush(stdout);
    }

    setUpSections();
  }

  // sections for cell and compute multicasts/reductions
  // iter -> cell -> send batch
  std::map<int, std::map<int, CProxySection_Comm> > sects;

  void setUpSections() {
    printf("section set up\n");
    //std::map<int, std::map<int, CommInfo> > commMap;
    for (std::map<int, std::map<int, CommInfo> >::iterator iter = commMap.begin();
         iter != commMap.end(); ++iter) {
      int iteration = iter->first;
      for (std::map<int, CommInfo>::iterator iter2 = iter->second.begin();
           iter2 != iter->second.end(); ++iter2) {
        int cellid = iter2->first;
        CommInfo& info = iter2->second;

        CkVec<CkArrayIndex1D> elms;
        for (std::list<SendTo>::iterator iter3 = info.sends.begin();
             iter3 != info.sends.end(); ++iter3) {
          elms.push_back(iter3->pe);
        }
        // @todo right now everyone creates and sets these up
        sects[iteration][cellid] =
          CProxySection_Comm::ckNew(commProxy.ckGetArrayID(),
                                    elms.getVec(),
                                    elms.size());
        //sects[iteration][cellid].warmup();
      }
    }
    //commMap[info.iter][linearize3D(sender)]
  }

  StateNode getNextState() {
  }
};

#endif /*__COMM_H__*/
