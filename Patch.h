#ifndef __PATCH_H__
#define __PATCH_H__

extern /* readonly */ CProxy_Main mainProxy;
extern /* readonly */ CProxy_Patch patchArray;
extern /* readonly */ CkGroupID mCastGrpID;
extern /* readonly */ double stepTime;
extern /* readonly */ int finalStepCount;

#include "ckmulticast.h"

//data message to be sent to computes
struct ParticleDataMsg : public CkMcastBaseMsg, public CMessage_ParticleDataMsg {
    vec3* part; //list of atoms
    int lengthAll;  //length of list
    int x;    //x coordinate of patch sending this message
    int y;    //y coordinate
    int z;    //z coordinate
    bool doAtSync;  //flag for indicating load balancing step

    //pack important information
    void pup(PUP::er &p){
    CMessage_ParticleDataMsg::pup(p);
      p | lengthAll;
      p | x; p | y; p | z;
      p | doAtSync;
      if (p.isUnpacking()){
        part = new vec3[lengthAll];
      }
      PUParray(p, part, lengthAll);
    } 
};

//chare used to represent a cell
class Patch : public CBase_Patch {
  private:
    Patch_SDAG_CODE   //SDAG code
      CkVec<Particle> particles;  //list of atoms
    int **computesList;   //my compute locations
    int stepCount;		// to count the number of steps, and decide when to stop
    int myNumParts;   //number of atoms in my cell
    bool done_lb;     //was load balancing done in last step?
    bool perform_lb;  //should I do load balancing in this step?
    int inbrs;        //number of interacting neighbors
    int updateCount;

    void migrateToPatch(Particle p, int &px, int &py, int &pz);
    void updateForce(double *forces, int lengthUp); //update forces after reduction
    void updateProperties();	//updates properties after receiving forces from computes
    void limitVelocity(Particle &p); //limit velcities to an upper limit
    Particle& wrapAround(Particle &p); //particles going out of right enters from left
    CProxySection_Compute mCastSecProxy; //handle to section proxy of computes
    void nextStep(); 

  public:
    Patch();
    Patch(CkMigrateMessage *msg);
    ~Patch();

    void createComputes();  //add my computes
    void createSection();   //created multicast section of computes
    void localCreateSection();
    void ftresume();
    void migrateParticles();
    void sendPositions();
    void ResumeFromSync();

    //pack important data when I move
    void pup(PUP::er &p){
      CBase_Patch::pup(p);
      __sdag_pup(p);
      p | particles;
      p | stepCount;		
      p | myNumParts;
      p | done_lb;
      p | perform_lb;
      p | updateCount;
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
      //adjust the multicast tree to give best performance after moving
      if (p.isUnpacking()){
        CkMulticastMgr *mg = CProxy_CkMulticastMgr(mCastGrpID).ckLocalBranch();
        mg->resetSection(mCastSecProxy);
        mg->setReductionClient(mCastSecProxy, new CkCallback(CkReductionTarget(Patch,reduceForces), thisProxy(thisIndex.x, thisIndex.y, thisIndex.z)));
      }
    } 

};

#endif
