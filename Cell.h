#ifndef __CELL_H__
#define __CELL_H__

#include "ckmulticast.h"
#include <string>

// Data message to be sent to computes
struct ParticleDataMsg : public CkMcastBaseMsg, public CMessage_ParticleDataMsg {
  vec3* part; // List of atoms
  int lengthAll; // Length of list
  int x, y, z; // (x,y,z) coordinate of cell sending this message

  ParticleDataMsg(const int x_, const int y_, const int z_, const int numPos)
    : x(x_), y(y_), z(z_), lengthAll(numPos) {}

  void pup(PUP::er &p) {
    CMessage_ParticleDataMsg::pup(p);
    p | lengthAll;
    p | x; p | y; p | z;
    if (p.isUnpacking()) part = new vec3[lengthAll];
    PUParray(p, part, lengthAll);
  }
};

class CellMap : public CkArrayMap {
private:
  int num_x, num_y, num_z, num_yz;
  float ratio;

public:
  CellMap(int num_x_, int num_y_, int num_z_) :
    num_x(num_x_), num_y(num_y_), num_z(num_z_) {
    num_yz = num_y * num_z;
    ratio = ((float)CkNumPes())/(num_x * num_yz);
  }

  CellMap(CkMigrateMessage *m) {}

  int registerArray(CkArrayIndex& numElements, CkArrayID aid) {
    return 0;
  }

  int procNum(int arrayHdl, const CkArrayIndex &idx) {
    int *index =(int *)idx.data();
    int patchID = index[2] + index[1] * num_z + index[0] * num_yz;
    return (int)(patchID * ratio);
  }
};

// Chare used to represent a cell
class Cell : public CBase_Cell {
private:
  Cell_SDAG_CODE;

  // List of atoms
  std::vector<Particle> particles;
  // My compute locations
  std::vector<CkArrayIndex6D> computesList;
  // To count the number of steps and decide when to stop
  int stepCount;
  // Number of atoms in my cell
  int myNumParts;
  // Number of interacting neighbors
  int inbrs;
  double stepTime;
  int updateCount;
  // Store kinetic energy - initial and final
  double energy[2];
  int numReadyCheckpoint;
  void migrateToCell(Particle p, int &px, int &py, int &pz);
  // Updates properties after receiving forces from computes
  void updateProperties(vec3 *forces);
  // Limit velocities to an upper limit
  void limitVelocity(Particle &p);
  // Particles going out of right enters from left
  Particle& wrapAround(Particle &p);
  // Handle to section proxy of computes
  CProxySection_Compute mCastSecProxy;

public:
  Cell();
  Cell(CkMigrateMessage *msg);
  ~Cell();
  void createComputes();
  void createSection();
  void migrateParticles();
  void sendPositions();
  void startCheckpoint(int);
  void pup(PUP::er &p);
};

#endif
