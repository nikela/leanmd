#ifndef __PHYSICS_H__
#define __PHYSICS_H__

#include "ckmulticast.h"

extern /* readonly */ CkGroupID mCastGrpID;

extern /* readonly */ int cellArrayDimX;	// Number of Chare Rows
extern /* readonly */ int cellArrayDimY;	// Number of Chare Columns
extern /* readonly */ int cellArrayDimZ;
extern /* readonly */ int finalStepCount; 
extern /* readonly */ CProxy_Cell cellArray;
extern /* readonly */ CProxy_Compute computeArray;
extern /*readonly*/ CProxy_Comm commProxy;

#define BLOCK_SIZE	512

//function to calculate forces among 2 lists of atoms
inline double calcPairForces(ParticleDataMsg* first, ParticleDataMsg* second,
                             int stepCount, CkIndex6D compute) {
  int i, j, jpart, ptpCutOffSqd, diff;
  int firstLen = first->lengthAll;
  int secondLen = second->lengthAll;
  double powTwenty, powTen, r, rsqd, f, fr;
  vec3 separation, force;
  double rSix, rTwelve;
  double energy = 0;
  int doEnergy = 0;
  if(stepCount == 1 || stepCount == finalStepCount)
    doEnergy = 1;

  vec3 *firstmsg = new vec3[firstLen+1];
  vec3 *secondmsg = new vec3[secondLen+1];
  //check for wrap around and adjust locations accordingly
  if (abs(first->x - second->x) > 1){
    diff = CELL_SIZE_X * cellArrayDimX;
    if (second->x < first->x)
      diff = -1 * diff; 
    for (i = 0; i < firstLen; i++)
      first->part[i].x += diff;
  }
  if (abs(first->y - second->y) > 1){
    diff = CELL_SIZE_Y * cellArrayDimY;
    if (second->y < first->y)
      diff = -1 * diff; 
    for (i = 0; i < firstLen; i++)
      first->part[i].y += diff;
  }
  if (abs(first->z - second->z) > 1){
    diff = CELL_SIZE_Z * cellArrayDimZ;
    if (second->z < first->z)
      diff = -1 * diff; 
    for (i = 0; i < firstLen; i++)
      first->part[i].z += diff;
  } 
  ptpCutOffSqd = PTP_CUT_OFF * PTP_CUT_OFF;
  powTen = pow(10.0, -10);
  powTwenty = pow(10.0, -20);

  int i1, j1;
  for(i1 = 0; i1 < firstLen; i1=i1+BLOCK_SIZE)
    for(j1 = 0; j1 < secondLen; j1=j1+BLOCK_SIZE)
      for(i = i1; i < i1+BLOCK_SIZE && i < firstLen; i++) {
        for(jpart = j1; jpart < j1+BLOCK_SIZE && jpart < secondLen; jpart++) {
          separation = first->part[i] - second->part[jpart];
          rsqd = dot(separation, separation);
          if (rsqd >= 0.001 && rsqd < ptpCutOffSqd) {
            rsqd = rsqd * powTwenty;
            r = sqrt(rsqd);
            rSix = ((double)rsqd) * rsqd * rsqd;
            rTwelve = rSix * rSix;
            f = (double)(VDW_A / rTwelve - VDW_B / rSix);
            if(doEnergy)
              energy += (double)( VDW_A / (12*rTwelve) - VDW_B / (6*rSix));
            fr = f /r;
            force = separation * (fr * powTen);
            firstmsg[i] += force;
            secondmsg[jpart] -= force;
          }
        }
      }

  int X = cellArrayDimX;
  int Y = cellArrayDimY;
  int Z = cellArrayDimZ;
  firstmsg[firstLen].x = first->z * Y * X + first->y * X + first->x;
  firstmsg[firstLen].y = 1;
  firstmsg[firstLen].z = 0;
  secondmsg[secondLen].x = second->z * Y * X + second->y * X + second->x;
  secondmsg[secondLen].y = 1;
  secondmsg[secondLen].z = 0;

  //CkMulticastMgr *mCastGrp = CProxy_CkMulticastMgr(mCastGrpID).ckLocalBranch();
  //CkGetSectionInfo(*mcast1, first);
  //mCastGrp->contribute(sizeof(vec3)*(firstLen+1), firstmsg, CkReduction::sum_double, *mcast1);
  //CkGetSectionInfo(*mcast2, second);
  //mCastGrp->contribute(sizeof(vec3)*(secondLen+1), secondmsg, CkReduction::sum_double, *mcast2);
  
  CkIndex3D t1;
  t1.x = first->x;
  t1.y = first->y;
  t1.z = first->z;
  commProxy[CkMyPe()].ckLocal()->depositForces(firstmsg, firstLen + 1, compute, stepCount, t1);

  CkIndex3D t2;
  t2.x = second->x;
  t2.y = second->y;
  t2.z = second->z;
  commProxy[CkMyPe()].ckLocal()->depositForces(secondmsg, secondLen + 1, compute, stepCount, t2);


  //delete [] firstmsg;
  //delete [] secondmsg;
  delete first;
  delete second;
  return energy;
}

//function to calculate forces among atoms in a single list
inline double calcInternalForces(ParticleDataMsg* first, int stepCount, CkIndex6D compute) {
  int i, j, ptpCutOffSqd;
  int firstLen = first->lengthAll;
  double powTwenty, powTen, firstx, firsty, firstz, rx, ry, rz, r, rsqd, fx, fy, fz, f, fr;
  vec3 firstpos, separation, force;
  double rSix, rTwelve;
  double energy = 0;
  int doEnergy = 0;
  if(stepCount == 1 || stepCount == finalStepCount)
    doEnergy = 1;
  vec3 *firstmsg = new vec3[firstLen + 1];

  ptpCutOffSqd = PTP_CUT_OFF * PTP_CUT_OFF;
  powTen = pow(10.0, -10);
  powTwenty = pow(10.0, -20);
  for(i = 0; i < firstLen; i++){
    firstpos = first->part[i];
    for(j = i+1; j < firstLen; j++) {
      // computing base values
      separation = firstpos - first->part[j];
      rsqd = dot(separation, separation);
      if(rsqd >= 0.001 && rsqd < ptpCutOffSqd){
        rsqd = rsqd * powTwenty;
        r = sqrt(rsqd);
        rSix = ((double)rsqd) * rsqd * rsqd;
        rTwelve = rSix * rSix;
        f = (double)(VDW_A / rTwelve - VDW_B / rSix);
        if(doEnergy)
          energy += (double)( VDW_A / (12*rTwelve) - VDW_B / (6*rSix));

        fr = f /r;
        force = separation * (fr * powTen);
        firstmsg[j] += force;
        firstmsg[i] -= force;
      }
    }
  }

  int X = cellArrayDimX;
  int Y = cellArrayDimY;
  firstmsg[firstLen].x = first->z * Y * X + first->y * X + first->x;
  firstmsg[firstLen].y = 1;
  firstmsg[firstLen].z = 0;

  //CkMulticastMgr *mCastGrp = CProxy_CkMulticastMgr(mCastGrpID).ckLocalBranch();
  //CkGetSectionInfo(*mcast1, first);
  //mCastGrp->contribute(sizeof(vec3)*(firstLen + 1), firstmsg, CkReduction::sum_double, *mcast1);

  CkIndex3D t1;
  t1.x = first->x;
  t1.y = first->y;
  t1.z = first->z;
  commProxy[CkMyPe()].ckLocal()->depositForces(firstmsg, firstLen + 1, compute, stepCount, t1);

  //delete [] firstmsg;
  delete first;
  return energy;
}
#endif
