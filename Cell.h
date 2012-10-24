#ifndef __CELL_H__
#define __CELL_H__
extern /* readonly */ CProxy_Main mainProxy;
extern /* readonly */ CProxy_Cell cellArray;
extern /* readonly */ CkGroupID mCastGrpID;
extern /* readonly */ int firstLdbStep;
extern /* readonly */ int ldbPeriod;
extern /* readonly */ int finalStepCount;
extern /* readonly */ int checkptFreq;

#include "ckmulticast.h"
#include <string>

//data message to be sent to computes
struct ParticleDataMsg : public CkMcastBaseMsg, public CMessage_ParticleDataMsg {
  vec3* part; //list of atoms
  int lengthAll;  //length of list
  int x, y, z;    // (x,y,z) coordinate of cell sending this message

  ParticleDataMsg(const int x_, const int y_, const int z_, const int numPos)
    : x(x_), y(y_), z(z_), lengthAll(numPos) { }

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
    std::vector<Particle> particles;  //list of atoms
    int **computesList;   //my compute locations
    int stepCount;		// to count the number of steps, and decide when to stop
    int myNumParts;   //number of atoms in my cell
    int inbrs;        //number of interacting neighbors
    double stepTime;
    int updateCount;
    double energy[2]; //store kinetic energy - initial and final
    int numReadyCheckpoint;
    void migrateToCell(Particle p, int &px, int &py, int &pz);
    void updateProperties(vec3 *forces);	//updates properties after receiving forces from computes
    void limitVelocity(Particle &p); //limit velcities to an upper limit
    Particle& wrapAround(Particle &p); //particles going out of right enters from left
    CProxySection_Compute mCastSecProxy; //handle to section proxy of computes

  public:
    Cell();
    Cell(CkMigrateMessage *msg);
    ~Cell();
    void createComputes();  //add my computes
    void createSection();   //created multicast section of computes
    void migrateParticles();
    void sendPositions();
    void startCheckpoint(int);
    void pup(PUP::er &p);
};

#endif
