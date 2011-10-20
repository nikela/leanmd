#ifndef __PHYSICS_H__
#define __PHYSICS_H__

#include "ckmulticast.h"

extern /* readonly */ CProxy_Main mainProxy;
extern /* readonly */ CProxy_Patch patchArray;
extern /* readonly */ CProxy_Compute computeArray;
extern /* readonly */ CkGroupID mCastGrpID;

extern /* readonly */ int patchArrayDimX;	// Number of Chare Rows
extern /* readonly */ int patchArrayDimY;	// Number of Chare Columns
extern /* readonly */ int patchArrayDimZ;
extern /* readonly */ int finalStepCount; 
extern /* readonly */ BigReal stepTime; 

#define BLOCK_SIZE	512

inline double calcPairForces(ParticleDataMsg* first, ParticleDataMsg* second, CkSectionInfo* cookie1, CkSectionInfo* cookie2, int stepCount) {
  int i, j, jpart, ptpCutOffSqd, diff;
  int firstLen = first->lengthAll;
  int secondLen = second->lengthAll;
  BigReal powTwenty, powTen, rx, ry, rz, r, rsqd, fx, fy, fz, f, fr;
  double rSix, rTwelve;
  double energy = 0;
  int doEnergy = 0;
  if(stepCount == 1 || stepCount == finalStepCount)
    doEnergy = 1;

  ParticleForceMsg *firstmsg = new (firstLen) ParticleForceMsg;
  ParticleForceMsg *secondmsg = new (secondLen) ParticleForceMsg;
  firstmsg->lengthUpdates = firstLen;
  secondmsg->lengthUpdates = secondLen;
  //check for wrap around and adjust locations accordingly
  if (abs(first->x - second->x) > 1){
    diff = PATCH_SIZE_X * patchArrayDimX;
    if (second->x < first->x)
      diff = -1 * diff; 
    for (i = 0; i < firstLen; i++)
      first->part[i].coord.x += diff;
  }
  if (abs(first->y - second->y) > 1){
    diff = PATCH_SIZE_Y * patchArrayDimY;
    if (second->y < first->y)
      diff = -1 * diff; 
    for (i = 0; i < firstLen; i++)
      first->part[i].coord.y += diff;
  }
  if (abs(first->z - second->z) > 1){
    diff = PATCH_SIZE_Z * patchArrayDimZ;
    if (second->z < first->z)
      diff = -1 * diff; 
    for (i = 0; i < firstLen; i++)
      first->part[i].coord.z += diff;
  } 
  ptpCutOffSqd = PTP_CUT_OFF * PTP_CUT_OFF;
  powTen = pow(10.0, -10);
  powTwenty = pow(10.0, -20);

  memset(firstmsg->forces, 0, firstLen * 3*sizeof(BigReal));
  memset(secondmsg->forces, 0, secondLen * 3*sizeof(BigReal));
  int i1, j1;
  for(i1 = 0; i1 < firstLen; i1=i1+BLOCK_SIZE)
    for(j1 = 0; j1 < secondLen; j1=j1+BLOCK_SIZE)
      for(i = i1; i < i1+BLOCK_SIZE && i < firstLen; i++) {
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
            f = (BigReal)(VDW_A / rTwelve - VDW_B / rSix);
            if(doEnergy)
	      energy += (BigReal)( VDW_A / (12*rTwelve) - VDW_B / (6*rSix));
            fr = f /r;
            fx = rx * fr * powTen;
            fy = ry * fr * powTen;
            fz = rz * fr * powTen;
            secondmsg->forces[jpart].x -= fx;
            secondmsg->forces[jpart].y -= fy;
            secondmsg->forces[jpart].z -= fz;
            firstmsg->forces[i].x += fx;
            firstmsg->forces[i].y += fy;
            firstmsg->forces[i].z += fz;
          }
        }
      }
  CkMulticastMgr *mCastGrp = CProxy_CkMulticastMgr(mCastGrpID).ckLocalBranch();
  CkGetSectionInfo(*cookie1, first);
  mCastGrp->contribute(sizeof(double)*3*firstmsg->lengthUpdates, firstmsg->forces, CkReduction::sum_double, *cookie1);
  CkGetSectionInfo(*cookie2, second);
  mCastGrp->contribute(sizeof(double)*3*secondmsg->lengthUpdates, secondmsg->forces, CkReduction::sum_double, *cookie2);

  delete firstmsg;
  delete secondmsg;
  delete first;
  delete second;
  return energy;
}

// Local function to compute all the interactions between pairs of particles in two sets
inline double calcInternalForces(ParticleDataMsg* first, CkSectionInfo *cookie1, int stepCount) {
  int i, j, ptpCutOffSqd;
  int firstLen = first->lengthAll;
  BigReal powTwenty, powTen, firstx, firsty, firstz, rx, ry, rz, r, rsqd, fx, fy, fz, f, fr;
  double rSix, rTwelve;
  double energy = 0;
  int doEnergy = 0;
  if(stepCount == 1 || stepCount == finalStepCount)
    doEnergy = 1;
  ParticleForceMsg *firstmsg = new (firstLen) ParticleForceMsg;
  firstmsg->lengthUpdates = firstLen;

  memset(firstmsg->forces, 1, firstLen * 3*sizeof(BigReal));
  ptpCutOffSqd = PTP_CUT_OFF * PTP_CUT_OFF;
  powTen = pow(10.0, -10);
  powTwenty = pow(10.0, -20);
  for(i = 0; i < firstLen; i++){
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
        f = (BigReal)(VDW_A / rTwelve - VDW_B / rSix);
        if(doEnergy)
	  energy += (BigReal)( VDW_A / (12*rTwelve) - VDW_B / (6*rSix));

        fr = f /r;
        fx = rx * fr * powTen;
        fy = ry * fr * powTen;
        fz = rz * fr * powTen;
        firstmsg->forces[j].x += fx;
        firstmsg->forces[j].y += fy;
        firstmsg->forces[j].z += fz;
        firstmsg->forces[i].x -= fx;
        firstmsg->forces[i].y -= fy;
        firstmsg->forces[i].z -= fz;
      }
    }
  }
  CkMulticastMgr *mCastGrp = CProxy_CkMulticastMgr(mCastGrpID).ckLocalBranch();
  mCastGrp->contribute(sizeof(double)*3*firstmsg->lengthUpdates, firstmsg->forces, CkReduction::sum_double, *cookie1);
  delete firstmsg;
  delete first;
  return energy;
}
#endif
