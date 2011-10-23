#include "time.h"
#include "defs.h"
#include "leanmd.decl.h"
#include "Patch.h"
#include "ckmulticast.h"

extern /* readonly */ CProxy_Main mainProxy;
extern /* readonly */ CProxy_Patch patchArray;
extern /* readonly */ CProxy_Compute computeArray;
extern /* readonly */ CkGroupID mCastGrpID;

extern /* readonly */ int patchArrayDimX;
extern /* readonly */ int patchArrayDimY;
extern /* readonly */ int patchArrayDimZ;
extern /* readonly */ int finalStepCount; 
extern /* readonly */ int firstLdbStep; 
extern /* readonly */ int ldbPeriod; 
extern /* readonly */ int ftPeriod; 
extern /* readonly */ double stepTime; 

//default constructor
Patch::Patch() {
  __sdag_init();
  int i;
  inbrs = NUM_NEIGHBORS;
  //total number of atoms in the system
  int numParts = PARTICLES_PER_PATCH * patchArrayDimX * patchArrayDimY * patchArrayDimZ;
  //load balancing to be called when AtSync is called
  usesAtSync = CmiTrue;

  myNumParts = PARTICLES_PER_PATCH;
  // starting random generator
  srand48(25763);

  // Particle initialization
  for(i=0; i < myNumParts; i++) {
    particles.push_back(Particle());
    particles[i].mass = HYDROGEN_MASS;

    //give random values for position and velocity
    particles[i].x = drand48() * PATCH_SIZE_X + thisIndex.x * PATCH_SIZE_X;
    particles[i].y = drand48() * PATCH_SIZE_Y + thisIndex.y * PATCH_SIZE_Y;
    particles[i].z = drand48() * PATCH_SIZE_Z + thisIndex.z * PATCH_SIZE_Z;
    particles[i].vx = (drand48() - 0.5) * .2 * MAX_VELOCITY;
    particles[i].vy = (drand48() - 0.5) * .2 * MAX_VELOCITY;
    particles[i].vz = (drand48() - 0.5) * .2 * MAX_VELOCITY;
    particles[i].fx = 0;
    particles[i].fy = 0;
    particles[i].fz = 0;

    particles[i].id = (thisIndex.x*patchArrayDimX + thisIndex.y) * numParts / (patchArrayDimX*patchArrayDimY)  + i;
  }

  stepCount = 0;
  done_lb = true;
  perform_lb = false;
}

//constructor for chare object migration
Patch::Patch(CkMigrateMessage *msg): CBase_Patch(msg) {
  __sdag_init();
  usesAtSync = CmiTrue;
  delete msg;
}  

Patch::~Patch() {}

//function to create my computes
void Patch::createComputes() {
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

  /*  The computes X are inserted by a given patch:
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
      px2 = x+dx+2;
      py1 = y + 2;
      py2 = y+dy+2;
      pz1 = z + 2;
      pz2 = z+dz+2;
      computeArray(px1, py1, pz1, px2, py2, pz2).insert((++currPe)%numPes);
      computesList[num][0] = px1; computesList[num][1] = py1; computesList[num][2] = pz1; 
      computesList[num][3] = px2; computesList[num][4] = py2; computesList[num][5] = pz2;
    }
    else {
      // these computes will be created by pairing patches
      px2 = WRAP_X(x+dx);
      py2 = WRAP_Y(y+dy);
      pz2 = WRAP_Z(z+dz);
      px1 = x;
      py1 = y;
      pz1 = z; 
      px1 = px2 - dx + 2;
      px2 = px2+2;
      py1 = py2 - dy + 2;
      py2 = py2+2;
      pz1 = pz2 - dz + 2;
      pz2 = pz2+2;
      computesList[num][0] = px2; computesList[num][1] = py2; computesList[num][2] = pz2; 
      computesList[num][3] = px1; computesList[num][4] = py1; computesList[num][5] = pz1;
    }
  } // end of for loop
  contribute(CkCallback(CkIndex_Main::startUpDone(), mainProxy));
}

//call multicast section creation
void Patch::createSection() {
  localCreateSection();
  contribute(CkCallback(CkIndex_Main::startUpDone(), mainProxy));
}

//function to create the multicast section of computes
void Patch::localCreateSection() {
  CkVec<CkArrayIndex6D> elems;
  //create a vector list of my computes
  for (int num=0; num<inbrs; num++)
    elems.push_back(CkArrayIndex6D(computesList[num][0], computesList[num][1], computesList[num][2], computesList[num][3], computesList[num][4], computesList[num][5]));

  CkArrayID computeArrayID = computeArray.ckGetArrayID();
  //knit the computes into a section
  mCastSecProxy = CProxySection_Compute::ckNew(computeArrayID, elems.getVec(), elems.size()); 

  //delegate the communication responsibility for this section to multicast library
  CkMulticastMgr *mCastGrp = CProxy_CkMulticastMgr(mCastGrpID).ckLocalBranch();
  mCastSecProxy.ckSectionDelegate(mCastGrp);
  mCastGrp->setReductionClient(mCastSecProxy, new CkCallback(CkReductionTarget(Patch,reduceForces), thisProxy(thisIndex.x, thisIndex.y, thisIndex.z)));
}

// Function to start interaction among particles in neighboring cells as well as its own particles
void Patch::sendPositions() {
  int len = particles.length();

  //create the particle and control message to be sent to computes
  ParticleDataMsg* msg = new (len) ParticleDataMsg;
  msg->x = thisIndex.x;
  msg->y = thisIndex.y;
  msg->z = thisIndex.z;
  msg->lengthAll = len;
  msg->doAtSync = false;

  //detect if load balancing is to be done
  if (stepCount >= firstLdbStep && (stepCount - firstLdbStep) % ldbPeriod == 0){
    msg->doAtSync = true;
    perform_lb = true;
  }
  for (int i = 0; i < len; i++){
    msg->part[i].x = particles[i].x;
    msg->part[i].y = particles[i].y;
    msg->part[i].z = particles[i].z;
  }
  mCastSecProxy.interact(msg);
}

//send the atoms that have moved beyond my cell to neighbors
void Patch::migrateParticles(){
  int i, x, y, z, x1, y1, z1;
  CkVec<Particle> *outgoing = new CkVec<Particle>[inbrs];

  // sending particles to neighboring cells
  x = thisIndex.x;
  y = thisIndex.y;
  z = thisIndex.z;

  for(i=0; i<particles.length(); i++) {
    migrateToPatch(particles[i], x1, y1, z1);
    if(x1 !=0 || y1!=0 || z1 !=0) {
      outgoing[(x1+1)*NBRS_Y*NBRS_Z + (y1+1)*NBRS_Z + (z1+1)].push_back(wrapAround(particles[i]));
      particles.remove(i);
    }
  }

  for(int num=0; num<inbrs; num++) {
    x1 = num / (NBRS_Y * NBRS_Z)            - NBRS_X/2;
    y1 = (num % (NBRS_Y * NBRS_Z)) / NBRS_Z - NBRS_Y/2;
    z1 = num % NBRS_Z                       - NBRS_Z/2;

    patchArray(WRAP_X(x+x1), WRAP_Y(y+y1), WRAP_Z(z+z1)).receiveParticles(outgoing[num]);
  }

  delete [] outgoing;
}

//check if the particle is to be moved
void Patch::migrateToPatch(Particle p, int &px, int &py, int &pz) {
  // currently this is assuming that particles are
  // migrating only to the immediate neighbors
  int x = thisIndex.x * PATCH_SIZE_X + PATCH_ORIGIN_X;
  int y = thisIndex.y * PATCH_SIZE_Y + PATCH_ORIGIN_Y;
  int z = thisIndex.z * PATCH_SIZE_Z + PATCH_ORIGIN_Z;

  if (p.x < x) px = -1;
  else if (p.x > x+PATCH_SIZE_X) px = 1;
  else px = 0;

  if (p.y < y) py = -1;
  else if (p.y > y+PATCH_SIZE_Y) py = 1;
  else py = 0;

  if (p.z < z) pz = -1;
  else if (p.z > z+PATCH_SIZE_Z) pz = 1;
  else pz = 0;
}

//function to decide on what to do now
void Patch::nextStep() {
  //if all steps are done, exit
  if(stepCount == finalStepCount) {
    contribute(CkCallback(CkIndex_Main::allDone(), mainProxy));
  } else if (perform_lb) {
    //call AtSync and go for load balancing
    perform_lb = false;
    done_lb = true;
    AtSync();
  } else {
    //continue to next step
    if(done_lb) {
      //reset timers if load balacing was just done
      done_lb = false;
      if((thisIndex.x + thisIndex.y + thisIndex.z) == 0)
        stepTime = CkWallTimer();
    }
    patchArray(thisIndex.x,thisIndex.y,thisIndex.z).doStep();
  }
}

//call nextStep if load balancing is done
void Patch::ResumeFromSync(){
  patchArray(thisIndex.x,thisIndex.y,thisIndex.z).nextStep();
}

//call nextStep if fault tolerance backup is done
  void Patch::ftresume(){
    if (thisIndex.x==0 && thisIndex.y==0 && thisIndex.z ==0)
      CkPrintf("patch 0 calling ftresume at %f\n",CmiWallTimer());
    patchArray(thisIndex.x,thisIndex.y,thisIndex.z).nextStep();
  }

//update forces on my atoms on values received from my computes
void Patch::updateForce(double *forces, int lengthUp) {
  int i;
  for(i = 0; i < lengthUp; i+=3) {
    particles[i/3].fx += forces[i];
    particles[i/3].fy += forces[i+1];
    particles[i/3].fz += forces[i+2];
  } 
}


// Function to update properties (i.e. acceleration, velocity and position) in particles
void Patch::updateProperties() {
  int i;
  double energy = 0;
  double powTen, powFteen, realTimeDelta, invMassParticle;
  powTen = pow(10.0, -10);
  powFteen = pow(10.0, -15);
  realTimeDelta = DEFAULT_DELTA * powFteen;
  for(i = 0; i < particles.length(); i++) {
    //calculate energy only in begining and end
    if(stepCount == 0 || stepCount == (finalStepCount - 1)) 
      energy += (0.5*particles[i].mass*(particles[i].vx*particles[i].vx + particles[i].vy*particles[i].vy*particles[i].vz*particles[i].vz));

    // applying kinetic equations
    invMassParticle = 1 / particles[i].mass;
    particles[i].ax = particles[i].fx * invMassParticle;
    particles[i].ay = particles[i].fy * invMassParticle;
    particles[i].az = particles[i].fz * invMassParticle;
    particles[i].vx = particles[i].vx + particles[i].ax * realTimeDelta;
    particles[i].vy = particles[i].vy + particles[i].ay * realTimeDelta;
    particles[i].vz = particles[i].vz + particles[i].az * realTimeDelta;

    limitVelocity(particles[i]);

    particles[i].x = particles[i].x + particles[i].vx * realTimeDelta;
    particles[i].y = particles[i].y + particles[i].vy * realTimeDelta;
    particles[i].z = particles[i].z + particles[i].vz * realTimeDelta;

    particles[i].fx = 0.0;
    particles[i].fy = 0.0;
    particles[i].fz = 0.0;
  }
  //reduction on energy only in begining and end
  if(stepCount == 0 || stepCount == (finalStepCount - 1)) 
    contribute(sizeof(double),&energy,CkReduction::sum_double,CkCallback(CkReductionTarget(Main,energySumK),mainProxy));
}

void Patch::limitVelocity(Particle &p) {
  if( fabs( p.vx ) > MAX_VELOCITY ) {
    if( p.vx < 0.0 )
      p.vx = -MAX_VELOCITY;
    else
      p.vx = MAX_VELOCITY;
  }

  if( fabs(p.vy) > MAX_VELOCITY ) {
    if( p.vy < 0.0 )
      p.vy = -MAX_VELOCITY;
    else
      p.vy = MAX_VELOCITY;
  }

  if( fabs(p.vz) > MAX_VELOCITY ) {
    if( p.vz < 0.0 )
      p.vz = -MAX_VELOCITY;
    else
      p.vz = MAX_VELOCITY;
  }
}

Particle& Patch::wrapAround(Particle &p) {
  if(p.x < PATCH_ORIGIN_X) p.x += PATCH_SIZE_X*patchArrayDimX;
  if(p.y < PATCH_ORIGIN_Y) p.y += PATCH_SIZE_Y*patchArrayDimY;
  if(p.z < PATCH_ORIGIN_Z) p.z += PATCH_SIZE_Z*patchArrayDimZ;
  if(p.x > PATCH_ORIGIN_X + PATCH_SIZE_X*patchArrayDimX) p.x -= PATCH_SIZE_X*patchArrayDimX;
  if(p.y > PATCH_ORIGIN_Y + PATCH_SIZE_Y*patchArrayDimY) p.y -= PATCH_SIZE_Y*patchArrayDimY;
  if(p.z > PATCH_ORIGIN_Z + PATCH_SIZE_Z*patchArrayDimZ) p.z -= PATCH_SIZE_Z*patchArrayDimZ;

  return p;
}


