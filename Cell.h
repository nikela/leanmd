#ifndef __CELL_H__
#define __CELL_H__

extern /* readonly */ CProxy_Main mainProxy;
extern /* readonly */ CProxy_Cell cellArray;
extern /* readonly */ CkGroupID mCastGrpID;
extern /* readonly */ int firstLdbStep;
extern /* readonly */ int ldbPeriod;
extern /* readonly */ int finalStepCount;

#include "ckmulticast.h"
#include "Main.h"

//data message to be sent to computes
struct ParticleDataMsg : public CkMcastBaseMsg, public CMessage_ParticleDataMsg {
  vec3* part; //list of atoms
  int lengthAll;  //length of list
  int x;    //x coordinate of cell sending this message
  int y;    //y coordinate
  int z;    //z coordinate

  void pup(PUP::er &p){
    CMessage_ParticleDataMsg::pup(p);
    p | lengthAll;
    p | x; p | y; p | z;
    if (p.isUnpacking()){
      part = new vec3[lengthAll];
    }
    PUParray(p, part, lengthAll);
  } 
};

//chare used to represent a cell
class Cell : public CBase_Cell {
  private:
    Cell_SDAG_CODE   //SDAG code
    CkVec<Particle> particles;  //list of atoms
    int **computesList;   //my compute locations
    int stepCount;		// to count the number of steps, and decide when to stop
    int myNumParts;   //number of atoms in my cell
    int inbrs;        //number of interacting neighbors
    int stepTime;
    int updateCount;
    double energy[2]; //store kinetic energy - initial and final

    // SDAG static scheduling state stuff
    int currentState;
    int stateCount;
    int thisCell;
    StateNode current;

    void migrateToCell(Particle p, int &px, int &py, int &pz);
    void updateProperties(vec3 *forces, int lengthUp);	//updates properties after receiving forces from computes
    void limitVelocity(Particle &p); //limit velcities to an upper limit
    Particle& wrapAround(Particle &p); //particles going out of right enters from left
    CProxySection_Compute mCastSecProxy; //handle to section proxy of computes

  public:
    Cell();
    Cell(CkMigrateMessage *msg);
    ~Cell();
    void pup(PUP::er &p);
    void createComputes();  //add my computes
    void createSection();   //created multicast section of computes
    void migrateParticles();
    void sendPositions();
};

#endif
