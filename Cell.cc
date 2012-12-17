#include "defs.h"
#include "leanmd.decl.h"
#include "Cell.h"
#include "ckmulticast.h"

Cell::Cell() : inbrs(NUM_NEIGHBORS), stepCount(1), updateCount(0), computesList(NUM_NEIGHBORS) {
  //load balancing to be called when AtSync is called
  usesAtSync = CmiTrue;

  int myid = thisIndex.z+cellArrayDimZ*(thisIndex.y+thisIndex.x*cellArrayDimY); 
  myNumParts = PARTICLES_PER_CELL_START + (myid*(PARTICLES_PER_CELL_END-PARTICLES_PER_CELL_START))/(cellArrayDimX*cellArrayDimY*cellArrayDimZ);

  // starting random generator
  srand48(myid);

  // Particle initialization
  for(int i = 0; i < myNumParts; i++) {
    particles.push_back(Particle());
    particles[i].mass = HYDROGEN_MASS;

    //uniformly place particles, avoid close distance among them
    particles[i].pos.x = (GAP/(float)2) + thisIndex.x * CELL_SIZE_X + ((i*KAWAY_Y*KAWAY_Z)/(PERDIM*PERDIM))*GAP;
    particles[i].pos.y = (GAP/(float)2) + thisIndex.y * CELL_SIZE_Y + (((i*KAWAY_Z)/PERDIM)%(PERDIM/KAWAY_Y))*GAP;
    particles[i].pos.z = (GAP/(float)2) + thisIndex.z * CELL_SIZE_Z + (i%(PERDIM/KAWAY_Z))*GAP; 
    //give random values for velocity
    particles[i].vel.x = (drand48() - 0.5) * .2 * MAX_VELOCITY;
    particles[i].vel.y = (drand48() - 0.5) * .2 * MAX_VELOCITY;
    particles[i].vel.z = (drand48() - 0.5) * .2 * MAX_VELOCITY;
  }
  energy[0] = energy[1] = 0;
  setMigratable(CmiFalse);
}

//constructor for chare object migration
Cell::Cell(CkMigrateMessage *msg): CBase_Cell(msg) {
  usesAtSync = CmiTrue;
  setMigratable(CmiFalse);
  delete msg;
}  

Cell::~Cell() {}

//function to create my computes
void Cell::createComputes() {
  int x = thisIndex.x, y = thisIndex.y, z = thisIndex.z;
  int px1, py1, pz1, dx, dy, dz, px2, py2, pz2;

  /*  The computes X are inserted by a given cell:
   *
   *	^  X  X  X
   *	|  0  X  X
   *	y  0  0  0
   *	   x ---->
   */

  // for round robin insertion
  int currPe = CkMyPe();

  for (int num = 0; num < inbrs; num++) {
    dx = num / (NBRS_Y * NBRS_Z)                - NBRS_X/2;
    dy = (num % (NBRS_Y * NBRS_Z)) / NBRS_Z     - NBRS_Y/2;
    dz = num % NBRS_Z                           - NBRS_Z/2;

    if (num >= inbrs / 2){
      px1 = x + KAWAY_X;
      py1 = y + KAWAY_Y;
      pz1 = z + KAWAY_Z;
      px2 = px1+dx;
      py2 = py1+dy;
      pz2 = pz1+dz;
      CkArrayIndex6D index(px1, py1, pz1, px2, py2, pz2);
      computeArray[index].insert((++currPe) % CkNumPes());
      computesList[num] = index;
    } else {
      // these computes will be created by pairing cells
      px1 = WRAP_X(x + dx) + KAWAY_X;
      py1 = WRAP_Y(y + dy) + KAWAY_Y;
      pz1 = WRAP_Z(z + dz) + KAWAY_Z;
      px2 = px1 - dx;
      py2 = py1 - dy;
      pz2 = pz1 - dz;
      CkArrayIndex6D index(px1, py1, pz1, px2, py2, pz2);
      computesList[num] = index;
    }
  } // end of for loop
  contribute(CkCallback(CkReductionTarget(Main,run),mainProxy));
}

//call multicast section creation
void Cell::createSection() {
  //knit the computes into a section
  mCastSecProxy = CProxySection_Compute::ckNew(computeArray.ckGetArrayID(), &computesList[0], computesList.size());

  //delegate the communication responsibility for this section to multicast library
  CkMulticastMgr *mCastGrp = CProxy_CkMulticastMgr(mCastGrpID).ckLocalBranch();
  mCastSecProxy.ckSectionDelegate(mCastGrp);
  mCastGrp->setReductionClient(mCastSecProxy, new CkCallback(CkReductionTarget(Cell,reduceForces), thisProxy(thisIndex.x, thisIndex.y, thisIndex.z)));
}

// Function to start interaction among particles in neighboring cells as well as its own particles
void Cell::sendPositions() {
  unsigned int len = particles.size();
  //create the particle and control message to be sent to computes
  ParticleDataMsg* msg = new (len) ParticleDataMsg(thisIndex.x, thisIndex.y, thisIndex.z, len);

  for(int i = 0; i < len; ++i)
    msg->part[i] = particles[i].pos;
  mCastSecProxy.calculateForces(msg);
}

//send the atoms that have moved beyond my cell to neighbors
void Cell::migrateParticles(){
  int x1, y1, z1;
  std::vector<std::vector<Particle> > outgoing;
  outgoing.resize(inbrs);

  for(std::vector<Particle>::iterator iter = particles.begin(); iter != particles.end(); ++iter) {
    migrateToCell(*iter, x1, y1, z1);
    if(x1!=0 || y1!=0 || z1!=0) {
      outgoing[(x1+KAWAY_X)*NBRS_Y*NBRS_Z + (y1+KAWAY_Y)*NBRS_Z + (z1+KAWAY_Z)].push_back(wrapAround(*iter));
      iter = particles.erase(iter);
    }
  }
  
  for(int num = 0; num < inbrs; num++) {
    x1 = num / (NBRS_Y * NBRS_Z)            - NBRS_X/2;
    y1 = (num % (NBRS_Y * NBRS_Z)) / NBRS_Z - NBRS_Y/2;
    z1 = num % NBRS_Z                       - NBRS_Z/2;
    cellArray(WRAP_X(thisIndex.x+x1), WRAP_Y(thisIndex.y+y1), WRAP_Z(thisIndex.z+z1)).receiveParticles(outgoing[num]);
  }
}

//check if the particle is to be moved
void Cell::migrateToCell(Particle p, int &px, int &py, int &pz) {
  int x = thisIndex.x * CELL_SIZE_X + CELL_ORIGIN_X;
  int y = thisIndex.y * CELL_SIZE_Y + CELL_ORIGIN_Y;
  int z = thisIndex.z * CELL_SIZE_Z + CELL_ORIGIN_Z;
  px = py = pz = 0;

  if (p.pos.x < (x-CELL_SIZE_X)) px = -2;
  else if (p.pos.x < x) px = -1;
  else if (p.pos.x > (x+2*CELL_SIZE_X)) px = 2;
  else if (p.pos.x > (x+CELL_SIZE_X)) px = 1;

  if (p.pos.y < (y-CELL_SIZE_Y)) py = -2;
  else if (p.pos.y < y) py = -1;
  else if (p.pos.y > (y+2*CELL_SIZE_Y)) py = 2;
  else if (p.pos.y > (y+CELL_SIZE_Y)) py = 1;

  if (p.pos.z < (z-CELL_SIZE_Z)) pz = -2;
  else if (p.pos.z < z) pz = -1;
  else if (p.pos.z > (z+2*CELL_SIZE_Z)) pz = 2;
  else if (p.pos.z > (z+CELL_SIZE_Z)) pz = 1;
}

// Function to update properties (i.e. acceleration, velocity and position) in particles
void Cell::updateProperties(vec3 *forces) {
  int i;
  double powTen, powTwenty, realTimeDeltaVel, invMassParticle;
  powTen = pow(10.0, 10);
  powTwenty = pow(10.0, -20);
  realTimeDeltaVel = DEFAULT_DELTA * powTwenty;
  for(i = 0; i < particles.size(); i++) {
    //calculate energy only in begining and end
    if(stepCount == 1) {
      energy[0] += (0.5 * particles[i].mass * dot(particles[i].vel, particles[i].vel) * powTen); // in milliJoules
    } else if(stepCount == finalStepCount) { 
      energy[1] += (0.5 * particles[i].mass * dot(particles[i].vel, particles[i].vel) * powTen);
    }
    // applying kinetic equations
    invMassParticle = 1 / particles[i].mass;
    particles[i].acc = forces[i] * invMassParticle; // in m/sec^2
    particles[i].vel += particles[i].acc * realTimeDeltaVel; // in A/fm

    limitVelocity(particles[i]);

    particles[i].pos += particles[i].vel * DEFAULT_DELTA; // in A
  }
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

void Cell::registerResumeClient(){
	CkCallback _cb(CkIndex_Cell::resumeFromChkp(),thisProxy(thisIndex.x,thisIndex.y,thisIndex.z));
	setChkpResumeClient(_cb);
}
//pack important data when I move/checkpoint
void Cell::pup(PUP::er &p) {
  
//  if(p.isChecking())
//  	CkPrintf("[%d][%d] pup cell\n",CmiMyPartition(),CkMyPe());
  
//  if(p.isChecking())
//  	p.skip();	  
  
  CBase_Cell::pup(p);
  __sdag_pup(p);
  
//  if(p.isChecking())
//  	CkPrintf("[%d][%d] pup particles\n",CmiMyPartition(),CkMyPe());
  p | particles;
  
//  if(p.isChecking())
//	  p.resume();
  
  p | stepCount;
//  if(p.isChecking())
//  	CkPrintf("[%d][%d] pup stepCount %d\n",CmiMyPartition(),CkMyPe(),stepCount);
  
//  if(p.isChecking())
//	  p.skip();
  p | myNumParts;
  p | updateCount;
  
//  if(p.isChecking())
//	  p.skipNext();
  p | stepTime;
  
  p | inbrs;
  p | numReadyCheckpoint;
  PUParray(p, energy, 2);

//  if(p.isChecking())
//    CkPrintf("[%d][%d] pup cell ends\n",CmiMyPartition(),CkMyPe());
  
  p | computesList;

  p | mCastSecProxy;
  
  //adjust the multicast tree to give best performance after moving
//  if (p.isUnpacking()||p.isChecking()){
  if (p.isUnpacking()){
    if(CkInRestarting()){
      createSection();
	  registerResumeClient();
    }
    else{
      CkMulticastMgr *mg = CProxy_CkMulticastMgr(mCastGrpID).ckLocalBranch();
      mg->resetSection(mCastSecProxy);
      mg->setReductionClient(mCastSecProxy, new CkCallback(CkReductionTarget(Cell,reduceForces), thisProxy(thisIndex.x, thisIndex.y, thisIndex.z)));
    }
  }
  
//  if(p.isChecking())
//	  p.skip();
}

