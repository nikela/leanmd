#ifndef __PATCH_H__
#define __PATCH_H__

extern /* readonly */ CProxy_Main mainProxy;
extern /* readonly */ CProxy_Patch patchArray;
extern /* readonly */ CkGroupID mCastGrpID;
extern /* readonly */ double stepTime;
extern /* readonly */ int finalStepCount;

#include "ckmulticast.h"

struct force {
  double x;
  double y;
  double z;
};

class loc{
  public:
    double x;
    double y;
    double z;

    void pup(PUP::er &p){
      p|x; p|y; p|z;
    }
};

class partData{
  public:
    loc coord;
    double charge;
    
    void pup(PUP::er &p){
      p|coord; p|charge; 
    }
};

class ParticleDataMsg : public CkMcastBaseMsg, public CMessage_ParticleDataMsg {
  public:
    partData* part;
    int lengthAll;
    int x;
    int y;
    int z;
    bool updateList;
    bool deleteList;
    bool doAtSync;

    void pup(PUP::er &p){
      CMessage_ParticleDataMsg::pup(p);
      p | lengthAll;
      p | x; p | y; p | z;
      p | updateList;
      p | deleteList;
      p | doAtSync;
      if (p.isUnpacking()){
	part = new partData[lengthAll];
      }
      PUParray(p, part, lengthAll);
    } 
};

/** \class Patch
 *  Class representing a cell in the grid. 
 *  We consider each cell as a square of LxL units
 */
class Patch : public CBase_Patch {
  private:
    Patch_SDAG_CODE;
    CkVec<Particle> particles;
    int **computesList;
    int forceCount;		// to count the returns from interactions
    int stepCount;		// to count the number of steps, and decide when to stop
    int updateCount;
    int myNumParts;
    bool done_lb;
    bool perform_lb;
    int inbrs;
    
 
    void migrateToPatch(Particle p, int &px, int &py, int &pz);
    void updateForce(double *forces, int lengthUp);
    void updateProperties();	// updates properties after receiving forces from computes
    void applyForces();
    void limitVelocity(Particle &p);
    Particle& wrapAround(Particle &p);
    CProxySection_Compute mCastSecProxy;
    void nextStep();

  public:
    Patch();
    Patch(CkMigrateMessage *msg);
    ~Patch();

    void createComputes();
    void createSection();
    void localCreateSection();
    void resume();
    void ftresume();
    void checkNextStep();	// checks whether to continue with next step
    void migrateParticles();
    void sendPositions();
    void ResumeFromSync();

    void pup(PUP::er &p){
      CBase_Patch::pup(p);
      __sdag_pup(p);
      p | particles;
      p | forceCount;		
      p | stepCount;		
      p | updateCount;
      p | myNumParts;
      p | done_lb;
      p | perform_lb;
      p | inbrs;

      if (p.isUnpacking()){
	computesList = new int*[inbrs];
	for (int i = 0; i < inbrs; i++){
	  computesList[i] = new int[6];
	}
      }

      for (int i = 0; i < inbrs; i++){
        PUParray(p, computesList[i], 6);
      }

      p | mCastSecProxy;
      if (p.isUnpacking()){
        CkMulticastMgr *mg = CProxy_CkMulticastMgr(mCastGrpID).ckLocalBranch();
        mg->resetSection(mCastSecProxy);
        mg->setReductionClient(mCastSecProxy, new CkCallback(CkReductionTarget(Patch,reduceForces), thisProxy(thisIndex.x, thisIndex.y, thisIndex.z)));
      }
    } 

};

#endif
