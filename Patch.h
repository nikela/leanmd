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

};

//chare used to represent a cell
class Patch : public CBase_Patch {
  private:
    Patch_SDAG_CODE;   //SDAG code
    CkVec<Particle> particles;  //list of atoms
    int **computesList;   //my compute locations
    int stepCount;		// to count the number of steps, and decide when to stop
    int myNumParts;   //number of atoms in my cell
    bool done_lb;     //was load balancing done in last step?
    bool perform_lb;  //should I do load balancing in this step?
    bool perform_ft;
    int inbrs;        //number of interacting neighbors
    int updateCount;
    int stepTime;

    void migrateToPatch(Particle p, int &px, int &py, int &pz);
    void updateProperties(vec3 *forces, int lengthUp);	//updates properties after receiving forces from computes
    void limitVelocity(Particle &p); //limit velcities to an upper limit
    Particle& wrapAround(Particle &p); //particles going out of right enters from left
    CProxySection_Compute mCastSecProxy; //handle to section proxy of computes
    void nextStep(); 

  public:
    Patch();
    Patch(CkMigrateMessage *msg);
    ~Patch();
    void pup(PUP::er &p);

    void createComputes();  //add my computes
    void createSection();   //created multicast section of computes
    void localCreateSection();
    void ftresume();
    void migrateParticles();
    void sendPositions();
    void ResumeFromSync();
};

#endif
