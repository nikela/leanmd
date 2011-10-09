#ifndef __PHYSICS_H__
#define __PHYSICS_H__

#include "pup.h"

typedef double BigReal;

#define AVAGADROS_NUMBER        (6.022141 * pow(10.0,23))
#define COULOMBS_CONSTANT       (8.987551 * pow(10.0,-9))
#define ELECTRON_CHARGE         (1.602176 * pow(10.0,-19))

#define NUM_PARTICLES		1000

#define USE_PAIRLISTS		true	// generally faster if true

#define DEFAULT_DELTA		1	// in femtoseconds

#define DEFAULT_FIRST_LDB	101
#define DEFAULT_LDB_PERIOD	500
#define DEFAULT_FT_PERIOD	100000

#define PATCHARRAY_DIM_X	3
#define PATCHARRAY_DIM_Y	3
#define PATCHARRAY_DIM_Z	3
#define PTP_CUT_OFF		13	//  cut off for atom to atom interactions
#define PATCH_MARGIN		0 	// constant difference between cut off and patch size
#define PATCH_SIZE_X		(PTP_CUT_OFF + PATCH_MARGIN)
#define PATCH_SIZE_Y		(PTP_CUT_OFF + PATCH_MARGIN)
#define PATCH_SIZE_Z		(PTP_CUT_OFF + PATCH_MARGIN)
#define PATCH_ORIGIN_X		0
#define PATCH_ORIGIN_Y		0
#define PATCH_ORIGIN_Z		0

#define MIGRATE_STEPCOUNT	5
#define DEFAULT_FINALSTEPCOUNT	101
#define MAX_VELOCITY		30.0

#define KAWAY_X			1
#define KAWAY_Y			1
#define KAWAY_Z			1
#define NBRS_X			(2*KAWAY_X+1)
#define NBRS_Y			(2*KAWAY_Y+1)
#define NBRS_Z			(2*KAWAY_Z+1)
#define NUM_NEIGHBORS		(NBRS_X * NBRS_Y * NBRS_Z)

#define WRAP_X(a)		(((a)+patchArrayDimX)%patchArrayDimX)
#define WRAP_Y(a)		(((a)+patchArrayDimY)%patchArrayDimY)
#define WRAP_Z(a)		(((a)+patchArrayDimZ)%patchArrayDimZ)

// Class for keeping track of the properties for a particle
class Particle {
  public:

    int id;
    BigReal mass;	// mass of the particle
    BigReal charge;     // charge of particle
    BigReal x;		// position in x axis
    BigReal y;		// position in y axis
    BigReal z;		// position in z axis

    BigReal fx;		// total forces on x axis
    BigReal fy;		// total forces on y axis
    BigReal fz;		// total forces on z axis

    BigReal ax;		// acceleration on x axis
    BigReal ay;		// acceleration on y axis
    BigReal az;		// acceleration on z axis

    BigReal vx;		// velocity on x axis
    BigReal vy;		// velocity on y axis
    BigReal vz;		// velocity on z axis

    // Default constructor
    Particle() {
      fx = fy = fz = 0.0;
    }

    // Function for pupping properties
    void pup(PUP::er &p) {
      p | id; p | mass; p | charge;
      p | x;  p | y;  p | z;
      p | fx; p | fy; p | fz;
      p | ax; p | ay; p | az;
      p | vx; p | vy; p | vz;
    }
};

extern /* readonly */ CProxy_Main mainProxy;
extern /* readonly */ CProxy_Patch patchArray;
extern /* readonly */ CProxy_Compute computeArray;
extern /* readonly */ CkGroupID mCastGrpID;

extern /* readonly */ bool usePairLists;
extern /* readonly */ int numParts;
extern /* readonly */ int patchArrayDimX;	// Number of Chare Rows
extern /* readonly */ int patchArrayDimY;	// Number of Chare Columns
extern /* readonly */ int patchArrayDimZ;
extern /* readonly */ int patchSizeX;
extern /* readonly */ int patchSizeY;
extern /* readonly */ int patchSizeZ;
extern /* readonly */ int ptpCutOff;
extern /* readonly */ int finalStepCount; 
extern /* readonly */ BigReal stepTime; 

extern /* readonly */ double A;		// Force Calculation parameter 1
extern /* readonly */ double B;		// Force Calculation parameter 2

#define BLOCK_SIZE	512

//calculate pair forces using pairlists
inline CkVec<int>* calcPairForcesPL(ParticleDataMsg* first, ParticleDataMsg* second, CkVec<int> *pairList, int *numLists, CkSectionInfo* cookie1, CkSectionInfo* cookie2) {
  int i, j, jpart, ptpCutOffSqd, diff;
  int firstLen = first->lengthAll;
  int secondLen = second->lengthAll;
  BigReal powTwenty, firstx, firsty, firstz, rx, ry, rz, r, rsqd, fx, fy, fz, f, fr, eField, constants;
  double rSix, rTwelve;
  ParticleForceMsg *firstmsg = new (firstLen) ParticleForceMsg;
  ParticleForceMsg *secondmsg = new (secondLen) ParticleForceMsg;
  firstmsg->lengthUpdates = firstLen;
  secondmsg->lengthUpdates = secondLen;
  //check for wrap around and adjust locations accordingly
  if (abs(first->x - second->x) > 1){
    diff = patchSizeX * patchArrayDimX;
    if (second->x < first->x)
      diff = -1 * diff; 
    for (i = 0; i < firstLen; i++)
      first->part[i].coord.x += diff;
  }
  if (abs(first->y - second->y) > 1){
    diff = patchSizeY * patchArrayDimY;
    if (second->y < first->y)
      diff = -1 * diff; 
    for (i = 0; i < firstLen; i++)
      first->part[i].coord.y += diff;
  }
  if (abs(first->z - second->z) > 1){
    diff = patchSizeZ * patchArrayDimZ;
    if (second->z < first->z)
      diff = -1 * diff; 
    for (i = 0; i < firstLen; i++)
      first->part[i].coord.z += diff;
  } 
  ptpCutOffSqd = ptpCutOff * ptpCutOff;
  powTwenty = pow(10.0, -20);

 //check if pairlist needs to be updated
  if (first->updateList){
    pairList = new CkVec<int>[firstLen];
    *numLists = firstLen;
    for (i = 0; i < firstLen; i++){
      firstx = first->part[i].coord.x;
      firsty = first->part[i].coord.y;
      firstz = first->part[i].coord.z;
      for (j = 0; j < secondLen; j++){
        rx = firstx - second->part[j].coord.x;
        ry = firsty - second->part[j].coord.y;
        rz = firstz - second->part[j].coord.z;
	rsqd = rx*rx + ry*ry + rz*rz;
	if (rsqd >= 0.001 && rsqd < ptpCutOffSqd)
	  pairList[i].push_back(j);
      }
    }
  }

  constants = COULOMBS_CONSTANT * ELECTRON_CHARGE * ELECTRON_CHARGE;

  int i1, j1, rootidx;
 
  memset(firstmsg->forces, 0, firstLen * 3*sizeof(BigReal));
  memset(secondmsg->forces, 0, secondLen * 3*sizeof(BigReal));
  for(i = 0; i < firstLen; i=i+BLOCK_SIZE)
   for(i1 = i; i1 < i+BLOCK_SIZE && i1 < firstLen; i1++) {
    eField = first->part[i1].charge;
    firstx = first->part[i1].coord.x;
    firsty = first->part[i1].coord.y;
    firstz = first->part[i1].coord.z;
    for(j = 0; j < pairList[i1].length(); j=j+BLOCK_SIZE)
     for(j1 = j; j1 < j+BLOCK_SIZE && j1 < pairList[i1].length(); j1++) {
      jpart = pairList[i1][j1];
      rx = firstx - second->part[jpart].coord.x;
      ry = firsty - second->part[jpart].coord.y;
      rz = firstz - second->part[jpart].coord.z;
      rsqd = rx*rx + ry*ry + rz*rz;
      if (rsqd >= 0.001){
	rsqd = rsqd * powTwenty;
	r = sqrt(rsqd);
	rSix = ((double)rsqd) * rsqd * rsqd;
	rTwelve = rSix * rSix;
        f = (BigReal)(A / rTwelve - B / rSix);
	f -= eField * constants * second->part[jpart].charge / rsqd;  //positive force should be attractive
	fr = f /r;
	fx = rx * fr;
	fy = ry * fr;
	fz = rz * fr;
	secondmsg->forces[jpart].x += fx;
	secondmsg->forces[jpart].y += fy;
	secondmsg->forces[jpart].z += fz;
	firstmsg->forces[i1].x -= fx;
	firstmsg->forces[i1].y -= fy;
	firstmsg->forces[i1].z -= fz;
      }
    }
  }
#ifdef USE_SECTION_MULTICAST
  CkMulticastMgr *mCastGrp = CProxy_CkMulticastMgr(mCastGrpID).ckLocalBranch();
  CkGetSectionInfo(*cookie1, first);
  mCastGrp->contribute(sizeof(BigReal)*3*firstmsg->lengthUpdates, firstmsg->forces, CkReduction::sum_double, *cookie1);
  CkGetSectionInfo(*cookie2, second);
  mCastGrp->contribute(sizeof(BigReal)*3*secondmsg->lengthUpdates, secondmsg->forces, CkReduction::sum_double, *cookie2);
  delete firstmsg;
  delete secondmsg;
#else
  patchArray(first->x, first->y, first->z).receiveForces(firstmsg);
  patchArray(second->x, second->y, second->z).receiveForces(secondmsg);
#endif
  if (first->deleteList){
    for (i = 0; i < firstLen; i++){
      pairList[i].removeAll();
    }
    delete [] pairList;
    *numLists = -1;
    pairList = NULL;
  }
  delete first;
  delete second;
  return pairList;
}

//calculate pair forces without using pairlists
inline void calcPairForces(ParticleDataMsg* first, ParticleDataMsg* second, CkSectionInfo* cookie1, CkSectionInfo* cookie2) {
  int i, j, jpart, ptpCutOffSqd, diff;
  int firstLen = first->lengthAll;
  int secondLen = second->lengthAll;
  BigReal powTwenty, rx, ry, rz, r, rsqd, fx, fy, fz, f, fr, eField, constants;
  double rSix, rTwelve;

  ParticleForceMsg *firstmsg = new (firstLen) ParticleForceMsg;
  ParticleForceMsg *secondmsg = new (secondLen) ParticleForceMsg;
  firstmsg->lengthUpdates = firstLen;
  secondmsg->lengthUpdates = secondLen;
  //check for wrap around and adjust locations accordingly
  if (abs(first->x - second->x) > 1){
    diff = patchSizeX * patchArrayDimX;
    if (second->x < first->x)
      diff = -1 * diff; 
    for (i = 0; i < firstLen; i++)
      first->part[i].coord.x += diff;
  }
  if (abs(first->y - second->y) > 1){
    diff = patchSizeY * patchArrayDimY;
    if (second->y < first->y)
      diff = -1 * diff; 
    for (i = 0; i < firstLen; i++)
      first->part[i].coord.y += diff;
  }
  if (abs(first->z - second->z) > 1){
    diff = patchSizeZ * patchArrayDimZ;
    if (second->z < first->z)
      diff = -1 * diff; 
    for (i = 0; i < firstLen; i++)
      first->part[i].coord.z += diff;
  } 
  ptpCutOffSqd = ptpCutOff * ptpCutOff;
  powTwenty = pow(10.0, -20);
  constants = COULOMBS_CONSTANT * ELECTRON_CHARGE * ELECTRON_CHARGE;
 
  memset(firstmsg->forces, 0, firstLen * 3*sizeof(BigReal));
  memset(secondmsg->forces, 0, secondLen * 3*sizeof(BigReal));
  int i1, j1;
  for(i1 = 0; i1 < firstLen; i1=i1+BLOCK_SIZE)
  for(j1 = 0; j1 < secondLen; j1=j1+BLOCK_SIZE)
   for(i = i1; i < i1+BLOCK_SIZE && i < firstLen; i++) {
    eField = first->part[i].charge;
    for(jpart = j1; jpart < j1+BLOCK_SIZE && jpart < secondLen; jpart++) {
      rx = first->part[i].coord.x - second->part[jpart].coord.x;
      ry = first->part[i].coord.y - second->part[jpart].coord.y;
      rz = first->part[i].coord.z - second->part[jpart].coord.z;
      rsqd = rx*rx + ry*ry + rz*rz;
      if (rsqd >= 0.001 && rsqd < ptpCutOffSqd){
	rsqd = rsqd * powTwenty;
	r = sqrt(rsqd);
	rSix = ((double)rsqd) * rsqd * rsqd;
	rTwelve = rSix * rSix;
        f = (BigReal)(A / rTwelve - B / rSix);
	f -= eField * constants * second->part[jpart].charge / rsqd;
	fr = f /r;
	fx = rx * fr;
	fy = ry * fr;
	fz = rz * fr;
	secondmsg->forces[jpart].x -= fx;
	secondmsg->forces[jpart].y -= fy;
	secondmsg->forces[jpart].z -= fz;
	firstmsg->forces[i].x += fx;
	firstmsg->forces[i].y += fy;
	firstmsg->forces[i].z += fz;
      }
    }
  }
#ifdef USE_SECTION_MULTICAST
  CkMulticastMgr *mCastGrp = CProxy_CkMulticastMgr(mCastGrpID).ckLocalBranch();
  CkGetSectionInfo(*cookie1, first);
  mCastGrp->contribute(sizeof(BigReal)*3*firstmsg->lengthUpdates, firstmsg->forces, CkReduction::sum_double, *cookie1);
  CkGetSectionInfo(*cookie2, second);
  mCastGrp->contribute(sizeof(BigReal)*3*secondmsg->lengthUpdates, secondmsg->forces, CkReduction::sum_double, *cookie2);
  delete firstmsg;
  delete secondmsg;
#else
  patchArray(first->x, first->y, first->z).receiveForces(firstmsg);
  patchArray(second->x, second->y, second->z).receiveForces(secondmsg);
#endif
  delete first;
  delete second;
}

// Local function to compute all the interactions between pairs of particles in two sets
inline void calcInternalForces(ParticleDataMsg* first, CkSectionInfo *cookie1) {
  int i, j, ptpCutOffSqd;
  int firstLen = first->lengthAll;
  BigReal powTwenty, firstx, firsty, firstz, rx, ry, rz, r, rsqd, fx, fy, fz, f, fr, eField, constants;
  double rSix, rTwelve;

  ParticleForceMsg *firstmsg = new (firstLen) ParticleForceMsg;
  firstmsg->lengthUpdates = firstLen;
    
  constants = COULOMBS_CONSTANT * ELECTRON_CHARGE * ELECTRON_CHARGE;
  
  memset(firstmsg->forces, 0, firstLen * 3*sizeof(BigReal));
  ptpCutOffSqd = ptpCutOff * ptpCutOff;
  powTwenty = pow(10.0, -20);
  for(i = 0; i < firstLen; i++){
    eField = first->part[i].charge;
    firstx = first->part[i].coord.x;
    firsty = first->part[i].coord.y;
    firstz = first->part[i].coord.z;
    for(j = i+1; j < firstLen; j++) {
      // computing base values
      rx = firstx - first->part[j].coord.x;
      ry = firsty - first->part[j].coord.y;
      rz = firstz - first->part[j].coord.z;
      rsqd = rx*rx + ry*ry + rz*rz;
      if(rsqd >= 0.001 && rsqd < ptpCutOffSqd){
	rsqd = rsqd * powTwenty;
	r = sqrt(rsqd);
	rSix = ((double)rsqd) * rsqd * rsqd;
	rTwelve = rSix * rSix;
        f = (BigReal)(A / rTwelve - B / rSix);
	f -= eField * constants * first->part[j].charge / (rsqd);
        
	fr = f /r;
	fx = rx * fr;
        fy = ry * fr;
        fz = rz * fr;
	firstmsg->forces[j].x += fx;
	firstmsg->forces[j].y += fy;
	firstmsg->forces[j].z += fz;
	firstmsg->forces[i].x -= fx;
	firstmsg->forces[i].y -= fy;
	firstmsg->forces[i].z -= fz;
      }
    }
  }
#ifdef USE_SECTION_MULTICAST
  CkMulticastMgr *mCastGrp = CProxy_CkMulticastMgr(mCastGrpID).ckLocalBranch();
  mCastGrp->contribute(sizeof(BigReal)*3*firstmsg->lengthUpdates, firstmsg->forces, CkReduction::sum_double, *cookie1);
  delete firstmsg;
#else
  patchArray(first->x, first->y, first->z).receiveForces(firstmsg);
#endif
  delete first;
}

