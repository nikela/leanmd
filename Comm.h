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
#include <vector>
#include <cassert>

extern /* readonly */ CProxy_Compute computeArray;

inline int linearize3D(CkIndex3D idx) {
  return idx.x * (cellArrayDimY + 1) * (cellArrayDimZ + 1) +
         idx.y * (cellArrayDimZ + 1) +
         idx.z;
}

// override for 6D and just pull necessary elements
inline int linearizeFake6D(CkIndex6D idx) {
  // just pull out values, repack and call the right method
  CkIndex3D nidx;
  nidx.x = idx.x1;
  nidx.y = idx.y1;
  nidx.z = idx.z1;
  return linearize3D(nidx);
}

inline CkIndex3D unlinearize3D(int cellid) {
  CkIndex3D indx;
  indx.x = cellid % (cellArrayDimZ + 1);
  indx.y = (cellid / (cellArrayDimZ + 1)) % (cellArrayDimY + 1);
  indx.z = (cellid / (cellArrayDimZ + 1) / (cellArrayDimY + 1)) % (cellArrayDimX + 1);
  return indx;
}

inline int linearize6D(CkIndex6D idx) {
  const int dim = 10;
  const int X = dim, Y = dim, Z = dim;
  const int X2 = dim, Y2 = dim, Z2 = dim;
  // @todo ??? fix this serialization
  return idx.x1 * Z2 * X  * X2 * Y  * Y2 +
         idx.x2 * X  * X2 * Y  * Y2 +
         idx.y1 * X2 * Y  * Y2 +
         idx.y2 * Y  * Y2 +
         idx.z1 * Y2 +
         idx.z2;
}

struct Startup : public CkMcastBaseMsg, public CMessage_Startup {
  int globalID;
  int iter;
  CkIndex6D compute;
  int reducingCell;
};

#include <cassert>

//# object types [ 0 -> CELL,
//#                1 -> COMPUTE,
//#                2 -> NONE ]

//# task types   [ 0 -> CELL_SEND,
//#                1 -> CELL_UPDATE,
//#                2 -> COMPUTE,
//#                3 -> MIGRATE ]

struct StaticSchedule : public CBase_StaticSchedule {
  std::map<int, std::vector<StateNode> > cellTransition;
  std::map<int, std::vector<StateNode> > computeTransition;

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
        case 4: node.thisObject.idx.y1 = x; break;
        case 5: node.thisObject.idx.z1 = x; break;
        case 6: node.thisObject.idx.x2 = x; break;
        case 7: node.thisObject.idx.y2 = x; break;
        case 8: node.thisObject.idx.z2 = x; break;
        case 9: node.taskType = x; break;
        case 10: node.pe = x; break;
        case 11: node.waitTask = x; break;
	case 12: node.relObject.objType = x; break;
        case 13: node.relObject.idx.x1 = x; break;
        case 14: node.relObject.idx.y1 = x; break;
        case 15: node.relObject.idx.z1 = x; break;
        case 16: node.relObject.idx.x2 = x; break;
        case 17: node.relObject.idx.y2 = x; break;
        case 18: node.relObject.idx.z2 = x;
          //printf("curVal = %d, objtype = %d\n", curVal, node.thisObject.objType);
          if (node.thisObject.objType == 0) {
            //printf("found cell\n");
            // cell
            CkIndex6D idx6 = node.thisObject.idx;
            CkIndex3D idx3; idx3.x = idx6.x1; idx3.y = idx6.y1; idx3.z = idx6.z1;
            cellTransition[linearize3D(idx3)].push_back(node);
          } else if (node.thisObject.objType == 1) {
            //printf("found compute\n");
            // compute
            computeTransition[linearize6D(node.thisObject.idx)].push_back(node);
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
        case 4: sender.y = x; break;
        case 5: sender.z = x; break;
        case 6: break;
        case 7: break;
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
                case 3: obj.idx.y1 = x; break;
                case 4: obj.idx.z1 = x; break;
                case 5: obj.idx.x2 = x; break;
                case 6: obj.idx.y2 = x; break;
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

  //commMap: iter <cell-obj> (spe,rpe) -> [ (pe,[obj,iter]) ]

  // sections for cell and compute multicasts/reductions
  // iter -> cell -> send batch
  std::map<int, std::map<int, CProxySection_Comm> > sects;

  void setUpSections() {
    // @todo at scale this might go to hell
    int globalID = CkMyPe() * 10000;
    printf("section set up\n");
    //std::map<int, std::map<int, CommInfo> > commMap;
    for (std::map<int, std::map<int, CommInfo> >::iterator iter = commMap.begin();
         iter != commMap.end(); ++iter) {
      int iteration = iter->first;
      for (std::map<int, CommInfo>::iterator iter2 = iter->second.begin();
           iter2 != iter->second.end(); ++iter2) {
        int cellid = iter2->first;
        CommInfo& info = iter2->second;

        if (info.rpe == CkMyPe() || info.spe == CkMyPe()) {
          CkVec<CkArrayIndex1D> elms;
          for (std::list<SendTo>::iterator iter3 = info.sends.begin();
               iter3 != info.sends.end(); ++iter3) {
            // @todo uncomment this once we are running on the right number of pes 
            elms.push_back(iter3->pe);
          }
          // @todo right now everyone creates and sets these up
          sects[iteration][cellid] =
            CProxySection_Comm::ckNew(commProxy.ckGetArrayID(),
                                      elms.getVec(),
                                      elms.size());
          CkMulticastMgr *mCastGrp = CProxy_CkMulticastMgr(mCastGrpID).ckLocalBranch();
          CkAssert(mCastGrp);
          sects[iteration][cellid].ckSectionDelegate(mCastGrp);

          Startup& start = *new Startup();
          if (info.rpe == CkMyPe())
            start.globalID = globalID;
          else
            start.globalID = -1;
          sects[iteration][cellid].warmup(&start);
        
          if (info.rpe == CkMyPe()) {
            mCastGrp->setReductionClient(sects[iteration][cellid],
                                         new CkCallback(CkReductionTarget(Comm,reduceForces),
                                                        commProxy[CkMyPe()]));

            for (std::list<SendTo>::iterator iter3 = info.sends.begin();
                 iter3 != info.sends.end(); ++iter3) {
              int pe = iter3->pe;

              //void warmupRed(int iter, int reducingCell, int globalID, int n, int* computes) {

              int iter = iteration;
              int reducingCell = cellid;
              int computes[iter3->sendTo.size()];
              int x = 0;

              for (std::list<Obj>::iterator iter4 = iter3->sendTo.begin();
                 iter4 != iter3->sendTo.end(); ++iter4) {
                computes[x] = linearize6D(iter4->idx);
                x++;
              }
              commProxy[pe].warmupRed(iter, reducingCell, globalID, iter3->sendTo.size(), computes);
            }
          }
          globalID++;
        }
      }
    }
  }

};

struct BufWarm {
  int iter;
  int* computes;
  int n;
  int reducingCell;
  int globalID;
};

class Comm : public CBase_Comm {
  private:
    std::set<int> myCells;
    std::set<int> myComputes;
    std::map< int, std::list< ParticleDataMsg* > > msgs;
    std::map< int, std::list< std::pair< vec3*, int> > > vecmsgs;
    CkIndex6D bufRelease;
    int bufReleaseType;

  public:
    Comm()
      : bufReleaseType(-1) { // at start there is not buffered release
    }

    Comm(CkMigrateMessage*) { }

    void registerCell(CkIndex3D indx) {
      int cellid = linearize3D(indx);
      //CkPrintf("%d: registering cell %d\n", CkMyPe(), cellid);
      myCells.insert(cellid);
      checkBufferedRelease();
      checkMsgs(indx);
    }

    void releaseNext(int type, CkIndex6D indx) {
      CkAssert(type == 0 || type == 1);
      if (type == 0) {
        // it's a cell
        int cellid = linearizeFake6D(indx);
        if (myCells.find(cellid) != myCells.end()) {
          CkPrintf("%d: releasing NEXT cell %d (%d,%d,%d), already on PE\n",
                   CkMyPe(), cellid, indx.x1, indx.y1, indx.z1);
          //CkAssert(bufReleaseType == -1 || bufRelease == indx);
          bufReleaseType = -1;
          // convert to 3d index
          CkIndex3D indx3;
          indx3.x = indx.x1;
          indx3.y = indx.y1;
          indx3.z = indx.z1;
          CkAssert(cellArray[indx3].ckLocal() != 0);
          cellArray[indx3].ckLocal()->commRelease();
        } else {
          CkPrintf("%d: trying to release, but cell %d not here yet\n",
                   CkMyPe(), cellid);
          bufferRelease(type, indx);
        }
      } else {
        assert(type == 1); // it's a compute
        int computeid = linearize6D(indx);
        if (myComputes.find(computeid) != myComputes.end()) {
          CkPrintf("%d: releasing NEXT compute %d, already on PE\n", CkMyPe(),
		   computeid);
          //CkAssert(bufReleaseType == -1 || bufRelease == indx);
          bufReleaseType = -1;
          computeArray[indx].ckLocal()->commRelease();
        } else {
          CkPrintf("%d: trying to release, but compute %d not here yet\n",
                   CkMyPe(), computeid);
          bufferRelease(type, indx);
        }
      }
    }

    void bufferRelease(int type, CkIndex6D indx) {
      // @protocol save this indx and wait for incoming
      // assert that there wasn't another one, or this is being recalled
      // become something arrived
      //CkAssert(bufReleaseType == -1 || indx == bufRelease);
      bufReleaseType = type;
      bufRelease = indx;
    }

    void checkBufferedRelease() {
      if (bufReleaseType != -1) releaseNext(bufReleaseType, bufRelease);
    }

    void registerCompute(CkIndex6D indx) {
      int computeid = linearize6D(indx);
      //CkPrintf("%d: registering compute %d\n", CkMyPe(), computeid);
      myComputes.insert(computeid);
      checkBufferedRelease();
      checkMsgs(indx);
    }

    void startMigrate(CkIndex3D indx, int pe) {
      cellArray[indx].ckLocal()->startMigrate(pe);
    }

    void startMigrateComp(CkIndex6D indx, int pe) {
      computeArray[indx].ckLocal()->startMigrate(pe);
    }

    void unregister(CkIndex3D indx) {
      int cellid = linearize3D(indx);
      assert(myCells.find(cellid) != myCells.end());
      myCells.erase(cellid);
    }

    void unregister(CkIndex6D indx) {
      int computeid = linearize6D(indx);
      assert(myComputes.find(computeid) != myComputes.end());
      myComputes.erase(computeid);
    }

    void tryDeliver(vec3* forces, int n) {
      int num = (int)forces[n-1].x;
      int count = (int)forces[n-1].y;
      // it's a sum, so divide the num by the number augmentations
      int cellid = num / count;
      CkIndex3D indx = unlinearize3D(cellid);
      if(myCells.find(num) != myCells.end()) {
        deliver(forces, n, indx);
      }
      else {
        vecmsgs[num].push_back(*(new std::pair<vec3*,int> (forces,n)));
      }
    }

    void tryDeliver(ParticleDataMsg* m) {
      /*CkPrintf("tryDeliver message = %p\n", m);*/
      CkIndex3D cindx;
      cindx.x = m->x;
      cindx.y = m->y;
      cindx.z = m->z;

      int iteration = m->iter;
      int cellid = linearize3D(cindx);
      bool found = false;
      std::list<CkIndex6D> indxs;

      StaticSchedule& sched = *staticSch.ckLocalBranch();

      for (std::list<SendTo>::iterator iter = sched.commMap[iteration][cellid].sends.begin();
           iter != sched.commMap[iteration][cellid].sends.end(); ++iter) {
        if (iter->pe == CkMyPe()) {
          found = true;
          for (std::list<Obj>::iterator iter2 = iter->sendTo.begin();
               iter2 != iter->sendTo.end(); ++iter2) {
            indxs.push_back(iter2->idx);
          }
        }
      }
      CkAssert(found);

      for (std::list<CkIndex6D>::iterator iter = indxs.begin(); iter != indxs.end(); ++iter) {
        int num = linearize6D(*iter);
        if (myComputes.find(num) != myComputes.end()) {
          CmiReference(UsrToEnv(m));
          deliver(m,*iter);
        } else {
          CmiReference(UsrToEnv(m));
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
      StaticSchedule& sched = *staticSch.ckLocalBranch();
      int cellid = linearize3D(cellIndx);
      CkAssert(cellArray[cellIndx].ckLocal() != 0);
      int iter = cellArray[cellIndx].ckLocal()->stepCount;
      CkPrintf("%d: sending positions: from cell %d at iter %d to some section\n",
	       CkMyPe(), cellid, iter);
      fflush(stdout);
      CkAssert(sched.sects[iter].find(cellid) != sched.sects[iter].end());
      sched.sects[iter][cellid].tryDeliver(m);
    }

    void checkMsgs(CkIndex3D indx) {
      int cellid = linearize3D(indx);
      for(std::list< std::pair<vec3*,int> >::iterator iter =
          vecmsgs[cellid].begin(); iter != vecmsgs[cellid].end(); ++iter) {
        deliver((*iter).first, (*iter).second, indx);
      }
      vecmsgs[cellid].clear();
    }

    void checkMsgs(CkIndex6D indx) {
      int num = linearize6D(indx);
      for(std::list< ParticleDataMsg* >::iterator iter = msgs[num].begin();
          iter != msgs[num].end(); ++iter) {
        deliver(*iter, indx);
      }
      msgs[num].clear();
    }

    std::map<int, CkSectionInfo> tempInfo;

    // iter -> cellid -> deposits -> section
    std::map<int, std::map<int, std::pair<int, CkSectionInfo> > > redInfo;
    // globalid -> iter -> computeid -> cellid
    std::map<int, BufWarm> tempInfo2;

    void warmup(Startup* s1) {
      if (s1->globalID != -1) {
        CkSectionInfo info;
        CkGetSectionInfo(info, s1);

        if (tempInfo2.find(s1->globalID) != tempInfo2.end()) {
          BufWarm b = tempInfo2[s1->globalID];
          CkAssert(redInfo[b.iter].find(b.reducingCell) == redInfo[b.iter].end());
          redInfo[b.iter][b.reducingCell] = std::make_pair(b.n, info);
          //   CkPrintf("%d: warmup extracting info, redNo = %d\n", info.get_redNo());
          //   fflush(stdout);
        } else {
          tempInfo[s1->globalID] = info;
        }
      }
    }

    void warmupRed(int iter, int reducingCell, int globalID, int n, int* computes) {
      if (tempInfo.find(globalID) == tempInfo.end()) {
        BufWarm buf;
        buf.iter = iter;
        buf.reducingCell = reducingCell;
        buf.globalID = globalID;
        buf.n = n;
        buf.computes = computes;
	tempInfo2[globalID] = buf;
      } else {
        CkAssert(redInfo[iter].find(reducingCell) == redInfo[iter].end());
        redInfo[iter][reducingCell] = std::make_pair(n, tempInfo[globalID]);
      }
    }

    // void warmupRed(Startup* sref) {
    //   Startup& s = *sref;
    //   if (tempInfo.find(s.globalID) == tempInfo.end()) {
    //     tempInfo2[s.globalID].push_back(sref);
    //   } else {
    //     CkAssert(tempInfo.find(s.globalID) != tempInfo.end());
    //     CkSectionInfo info = tempInfo[s.globalID];
    //     CkAssert(redInfo[s.iter][linearize6D(s.compute)].find(s.reducingCell) ==
    //     	 redInfo[s.iter][linearize6D(s.compute)].end());
    //     redInfo[s.iter][linearize6D(s.compute)][s.reducingCell] = info;
    //     delete &s;
    //   }
    // }

    std::map<int, std::map<int, int> > redCount;
    std::map<int, std::map<int, std::pair<int, vec3*> > > redCombine;

    void depositForces(vec3* msg, int n, CkIndex6D depositor, int iter, CkIndex3D target) {
      int computeid = linearize6D(depositor);
      int cellid = linearize3D(target);
      fflush(stdout);
      CkAssert(redInfo[iter].find(cellid) != redInfo[iter].end());
      std::pair<int, CkSectionInfo>& p = redInfo[iter][cellid];
      CkSectionInfo info = p.second;

      CkPrintf("%d: depositing positions: from compute %d to cell %d at iter %d to"
               " section redNo = %d, val = %p\n",
	       CkMyPe(), computeid, linearize3D(target), iter, info.get_redNo(), info.get_val());
      fflush(stdout);

      if (redCombine[iter].find(cellid) == redCombine[iter].end()) {
        redCombine[iter][cellid] = std::make_pair(n, msg);
      } else {
        std::pair<int, vec3*>& payload = redCombine[iter][cellid];
        CkAssert(payload.first == n);
        for (int i = 0; i < n; i++) {
          vec3& vals = *payload.second;
          vals += *msg;
        }
      }
      
      redCount[iter][cellid]++;

      if (redCount[iter][cellid] == p.first) {
        CkMulticastMgr *mCastGrp = CProxy_CkMulticastMgr(mCastGrpID).ckLocalBranch();
        CkAssert(redCombine[iter][cellid].first == n);
        mCastGrp->contribute(sizeof(vec3)*n, redCombine[iter][cellid].second,
                             CkReduction::sum_double, info);
        CkPrintf("%d: contributing to reduction on elements redNo = %d, val = %p\n",
                 CkMyPe(), info.get_redNo(), info.get_val());
      }
    }
};

#endif /*__COMM_H__*/
