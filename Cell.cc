
#include "defs.h"
#include "leanmd.decl.h"
#include "Cell.h"
#include "Comm.h"
#include "ckmulticast.h"

extern /* readonly */ CProxy_PME pmeArray;
extern /* readonly */ CProxy_Main mainProxy;
extern /* readonly */ CProxy_Cell cellArray;
extern /* readonly */ CProxy_Compute computeArray;
extern /* readonly */ CkGroupID mCastGrpID;
extern /* readonly */ int cellArrayDimX;
extern /* readonly */ int cellArrayDimY;
extern /* readonly */ int cellArrayDimZ;
extern /* readonly */ int finalStepCount; 

//default constructor
Cell::Cell() {
  __sdag_init();
  int i;
  inbrs = NUM_NEIGHBORS;
  //load balancing to be called when AtSync is called
  usesAtSync = CmiTrue;

  myNumParts = PARTICLES_PER_CELL;
  // starting random generator
  srand48(thisIndex.x+cellArrayDimX*(thisIndex.y+thisIndex.z*cellArrayDimY));

  // Particle initialization
  for(i=0; i < myNumParts; i++) {
    particles.push_back(Particle());
    particles[i].mass = HYDROGEN_MASS;

    //give random values for position and velocity
    particles[i].pos.x = drand48() * CELL_SIZE_X + thisIndex.x * CELL_SIZE_X;
    particles[i].pos.y = drand48() * CELL_SIZE_Y + thisIndex.y * CELL_SIZE_Y;
    particles[i].pos.z = drand48() * CELL_SIZE_Z + thisIndex.z * CELL_SIZE_Z;
    particles[i].vel.x = (drand48() - 0.5) * .2 * MAX_VELOCITY;
    particles[i].vel.y = (drand48() - 0.5) * .2 * MAX_VELOCITY;
    particles[i].vel.z = (drand48() - 0.5) * .2 * MAX_VELOCITY;
  }

  stepCount = 1;
  updateCount = 0;
  stepTime = 0; 
  energy[0] = energy[1] = 0;
  
  commProxy[CkMyPe()].ckLocal()->registerCell(thisIndex);
}

//constructor for chare object migration
Cell::Cell(CkMigrateMessage *msg): CBase_Cell(msg) {
  __sdag_init();
  usesAtSync = CmiTrue;
  delete msg;
}  

Cell::~Cell() {}

//function to create my computes
void Cell::createComputes() {
  int num;  

  int x = thisIndex.x;
  int y = thisIndex.y;
  int z = thisIndex.z;
  int px1, py1, pz1, dx, dy, dz, px2, py2, pz2;

  // For Round Robin insertion
  int numPes = CkNumPes();
  int currPe = CkMyPe();

  computesList = new int*[inbrs];
  for (int i =0; i < inbrs; i++){
    computesList[i] = new int[6];
  }

  /*  The computes X are inserted by a given cell:
   *
   *	^  X  X  X
   *	|  0  X  X
   *	y  0  0  0
   *	   x ---->
   */

  for (num=0; num<inbrs; num++) {
    dx = num / (NBRS_Y * NBRS_Z)            - NBRS_X/2;
    dy = (num % (NBRS_Y * NBRS_Z)) / NBRS_Z   - NBRS_Y/2;
    dz = num % NBRS_Z                       - NBRS_Z/2;

    if (num >= inbrs/2){
      px1 = x + 2;
      py1 = y + 2;
      pz1 = z + 2;
      px2 = px1+dx;
      py2 = py1+dy;
      pz2 = pz1+dz;
      computeArray(px1, py1, pz1, px2, py2, pz2).insert();
      computesList[num][0] = px1; computesList[num][1] = py1; computesList[num][2] = pz1; 
      computesList[num][3] = px2; computesList[num][4] = py2; computesList[num][5] = pz2;
    }
    else {
      // these computes will be created by pairing cells
      px1 = WRAP_X(x+dx)+2 ;
      py1 = WRAP_Y(y+dy)+2;
      pz1 = WRAP_Z(z+dz)+2;
      px2 = px1 - dx;
      py2 = py1 - dy;
      pz2 = pz1 - dz;
      computesList[num][0] = px1; computesList[num][1] = py1; computesList[num][2] = pz1; 
      computesList[num][3] = px2; computesList[num][4] = py2; computesList[num][5] = pz2;
    }
  } // end of for loop
  contribute(0,NULL,CkReduction::nop,CkCallback(CkReductionTarget(Main,computesCreated),mainProxy));
}

//call multicast section creation
void Cell::createSection() {
  CkVec<CkArrayIndex6D> elems;
  //create a vector list of my computes
  for (int num=0; num<inbrs; num++) {
    elems.push_back(CkArrayIndex6D(computesList[num][0], computesList[num][1], computesList[num][2], computesList[num][3], computesList[num][4], computesList[num][5]));
    fflush(stdout);
  }

  CkArrayID computeArrayID = computeArray.ckGetArrayID();
  //knit the computes into a section
  mCastSecProxy = CProxySection_Compute::ckNew(computeArrayID, elems.getVec(), elems.size()); 

  //delegate the communication responsibility for this section to multicast library
  CkMulticastMgr *mCastGrp = CProxy_CkMulticastMgr(mCastGrpID).ckLocalBranch();
  mCastSecProxy.ckSectionDelegate(mCastGrp);
  mCastGrp->setReductionClient(mCastSecProxy, new CkCallback(CkReductionTarget(Comm,reduceForces), thisProxy(thisIndex.x, thisIndex.y, thisIndex.z)));
}

// Function to start interaction among particles in neighboring cells as well as its own particles
void Cell::sendPositions() {
  int len = particles.length();

  //create the particle and control message to be sent to computes
  ParticleDataMsg* msg = new (len) ParticleDataMsg;
  msg->x = thisIndex.x;
  msg->y = thisIndex.y;
  msg->z = thisIndex.z;
  msg->lengthAll = len;
  msg->iter = stepCount;

  for (int i = 0; i < len; i++)
    msg->part[i] = particles[i].pos;

  //PREMPTING WITH DAG SCHEDULING
  //mCastSecProxy.calculateForces(msg);
  commProxy[CkMyPe()].ckLocal()->calculateForces(msg);
}

//send the atoms that have moved beyond my cell to neighbors
void Cell::migrateParticles(){
}

//check if the particle is to be moved
void Cell::migrateToCell(Particle p, int &px, int &py, int &pz) {
  // currently this is assuming that particles are
  // migrating only to the immediate neighbors
  int x = thisIndex.x * CELL_SIZE_X + CELL_ORIGIN_X;
  int y = thisIndex.y * CELL_SIZE_Y + CELL_ORIGIN_Y;
  int z = thisIndex.z * CELL_SIZE_Z + CELL_ORIGIN_Z;
  px = py = pz = 0;

  if (p.pos.x < x) px = -1;
  else if (p.pos.x > x+CELL_SIZE_X) px = 1;

  if (p.pos.y < y) py = -1;
  else if (p.pos.y > y+CELL_SIZE_Y) py = 1;

  if (p.pos.z < z) pz = -1;
  else if (p.pos.z > z+CELL_SIZE_Z) pz = 1;
}

// Function to update properties (i.e. acceleration, velocity and position) in particles
void Cell::updateProperties(vec3 *forces, int lengthUp) {
  int i;
  double powTen, powFteen, realTimeDelta, invMassParticle;
  powTen = pow(10.0, -10);
  powFteen = pow(10.0, -15);
  realTimeDelta = DEFAULT_DELTA * powFteen;
  for(i = 0; i < particles.length(); i++) {
    //calculate energy only in begining and end
    if(stepCount == 1) {
      energy[0] += (0.5 * particles[i].mass * dot(particles[i].vel, particles[i].vel));
    } else if(stepCount == finalStepCount) { 
      energy[1] += (0.5 * particles[i].mass * dot(particles[i].vel, particles[i].vel));
    }
    // applying kinetic equations
    invMassParticle = 1 / particles[i].mass;
    particles[i].acc = forces[i] * invMassParticle;
    particles[i].vel += particles[i].acc * realTimeDelta;

    limitVelocity(particles[i]);

    particles[i].pos += particles[i].vel * realTimeDelta;
  }
}

void Cell::sendCharges() {
  double charges[CHARGES_PER_CELL];
  commProxy[CkMyPe()].ckLocal()->sendCharges(CHARGES_PER_CELL, charges,
					     thisIndex.x, thisIndex.y, stepCount, 0);
  //pmeArray(thisIndex.x,thisIndex.y).recvCharges(3000, charges);
}

inline double velocityCheck(double inVelocity) {
  if(fabs(inVelocity) > MAX_VELOCITY) {
    if(inVelocity < 0.0 )
      return -MAX_VELOCITY;
    else
      return MAX_VELOCITY;
  } else {
    return inVelocity;
  }
}

void Cell::limitVelocity(Particle &p) {
  p.vel.x = velocityCheck(p.vel.x);
  p.vel.y = velocityCheck(p.vel.y);
  p.vel.z = velocityCheck(p.vel.z);
}

Particle& Cell::wrapAround(Particle &p) {
  if(p.pos.x < CELL_ORIGIN_X) p.pos.x += CELL_SIZE_X*cellArrayDimX;
  if(p.pos.y < CELL_ORIGIN_Y) p.pos.y += CELL_SIZE_Y*cellArrayDimY;
  if(p.pos.z < CELL_ORIGIN_Z) p.pos.z += CELL_SIZE_Z*cellArrayDimZ;
  if(p.pos.x > CELL_ORIGIN_X + CELL_SIZE_X*cellArrayDimX) p.pos.x -= CELL_SIZE_X*cellArrayDimX;
  if(p.pos.y > CELL_ORIGIN_Y + CELL_SIZE_Y*cellArrayDimY) p.pos.y -= CELL_SIZE_Y*cellArrayDimY;
  if(p.pos.z > CELL_ORIGIN_Z + CELL_SIZE_Z*cellArrayDimZ) p.pos.z -= CELL_SIZE_Z*cellArrayDimZ;

  return p;
}

//pack important data when I move
void Cell::pup(PUP::er &p) {
  CBase_Cell::pup(p);
  __sdag_pup(p);
  p | particles;
  p | stepCount;
  p | myNumParts;
  p | updateCount;
  p | stepTime;
  p | inbrs;
  p | currentState;
  p | stateCount;
  p | thisCell;
  p | current;
  PUParray(p, energy, 2);

  if (p.isUnpacking()){
    computesList = new int*[inbrs];
    for (int i = 0; i < inbrs; i++){
      computesList[i] = new int[6];
    }
  }

  for (int i = 0; i < inbrs; i++){
    PUParray(p, computesList[i], 6);
  }
}

void Cell::startMigrate(int pe) {
  //CkPrintf("%d: cell startMigration to %d\n", CkMyPe(), pe);
  migrateMe(pe);
}

void Cell::ckJustMigrated() {
  //CkPrintf("%d: cell just ckJustMigrated on %d\n", CkMyPe(), thisCell);
  //fflush(stdout);
  ArrayElement::ckJustMigrated();
  migrateDone(0);
}


PME::PME()
  : charges(CHARGES_PER_CELL * cellArrayDimZ) {
  __sdag_init();
  
  CkIndex3D idx3;
  idx3.x = thisIndex.x;
  idx3.y = thisIndex.y;
  idx3.z = 0;
  commProxy[CkMyPe()].ckLocal()->registerPME(idx3);
}
