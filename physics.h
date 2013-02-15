#ifndef __PHYSICS_H__
#define __PHYSICS_H__

#include <xmmintrin.h>
#include <smmintrin.h>

extern /* readonly */ CkGroupID mCastGrpID;

extern /* readonly */ int cellArrayDimX;	// Number of Chare Rows
extern /* readonly */ int cellArrayDimY;	// Number of Chare Columns
extern /* readonly */ int cellArrayDimZ;
extern /* readonly */ int finalStepCount; 

#define BLOCK_SIZE	512

//function to calculate forces among 2 lists of atoms
inline double calcPairForces(ParticleDataMsg* first, ParticleDataMsg* second, int stepCount) {
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

  vec3 *firstmsg = new vec3[firstLen];
  vec3 *secondmsg = new vec3[secondLen];
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

//sse begin

double __attribute__((aligned(16))) first_x[firstLen], first_y[firstLen], first_z[firstLen];
double __attribute__((aligned(16))) sec_x[firstLen], sec_y[firstLen], sec_z[firstLen];
double __attribute__((aligned(16))) firstmsg_x[firstLen], firstmsg_y[firstLen], firstmsg_z[firstLen];
double __attribute__((aligned(16))) secmsg_x[firstLen], secmsg_y[firstLen], secmsg_z[firstLen];

//copy the values to the temp arrays
for(int i=0; i<firstLen; i++){
	first_x[i] = (first->part[i]).x;
	first_y[i] = (first->part[i]).y;
	first_z[i] = (first->part[i]).z;
	firstmsg_x[i] = firstmsg[i].x;
	firstmsg_y[i] = firstmsg[i].y;
	firstmsg_z[i] = firstmsg[i].z;

}
for(int i=0; i<secondLen; i++){
	sec_x[i] = (second->part[i]).x;
	sec_y[i] = (second->part[i]).y;
	sec_z[i] = (second->part[i]).z;
	secmsg_x[i] = secondmsg[i].x;
	secmsg_y[i] = secondmsg[i].y;
	secmsg_z[i] = secondmsg[i].z;
}
    
__m128d f_x, f_y, f_z, s_x, s_y, s_z;
__m128d fm_x, fm_y, fm_z, sm_x, sm_y, sm_z;
__m128d sp_x, sp_y, sp_z, rsqd_x, rsqd_y, rsqd_z;
__m128d force_x, force_y, force_z;

  int i1, j1;
//ckout<<": "<<firstLen<<endl;
  for(i1 = 0; i1 < firstLen; i1=i1+BLOCK_SIZE)
    for(j1 = 0; j1 < secondLen; j1=j1+BLOCK_SIZE)
      for(i = i1; i < i1+BLOCK_SIZE && i < firstLen; i+=2) {
        for(jpart = j1; jpart < j1+BLOCK_SIZE && jpart < secondLen; jpart+=2) {
            //ckout<<": "<<i1<<endl;
        //load the arrays to vectors
		f_x = _mm_load_pd(&first_x[i]);
		f_y = _mm_load_pd(&first_y[i]);
		f_z = _mm_load_pd(&first_z[i]);
		s_x = _mm_load_pd(&sec_x[jpart]);
		s_y = _mm_load_pd(&sec_y[jpart]);
		s_z = _mm_load_pd(&sec_z[jpart]);
        fm_x = _mm_load_pd(&firstmsg_x[i]);
        fm_y = _mm_load_pd(&firstmsg_y[i]);
        fm_z = _mm_load_pd(&firstmsg_z[i]);
        sm_x = _mm_load_pd(&secmsg_x[jpart]);
        sm_y = _mm_load_pd(&secmsg_y[jpart]);
        sm_z = _mm_load_pd(&secmsg_z[jpart]);
            
        sp_x = _mm_sub_pd(f_x, s_x); sp_y = _mm_sub_pd(f_y, s_y); sp_z = _mm_sub_pd(f_z, s_z);
        rsqd_x = _mm_mul_pd(sp_x, sp_x); rsqd_y = _mm_mul_pd(sp_y, sp_y); rsqd_z = _mm_mul_pd(sp_z, sp_z); //dot
        
        double rsqd0, rsqd1, rsqd2, rsqd3;
        rsqd0 = _mm_cvtsd_f64(rsqd_x) + _mm_cvtsd_f64(rsqd_y) + _mm_cvtsd_f64(rsqd_z);
        rsqd1 = _mm_cvtsd_f64(rsqd_x) + _mm_cvtsd_f64(rsqd_y) + _mm_cvtsd_f64(rsqd_z);
        
        force_x = _mm_set1_pd(0.0); force_y = _mm_set1_pd(0.0); force_z = _mm_set1_pd(0.0);
        
        if (rsqd0 >= 0.001 && rsqd0 < ptpCutOffSqd) {
            rsqd0 = rsqd0 * powTwenty;
            r = sqrt(rsqd0);
            rTwelve = rSix * rSix;
            rSix = ((double)rsqd0) * rsqd0 * rsqd0;
            f = (double)(VDW_A / rTwelve - VDW_B / rSix);
			if(doEnergy){
                energy += (double)( VDW_A / (12*rTwelve) - VDW_B / (6*rSix));
			}
            fr = f /r;

            force_x = _mm_castsi128_pd (_mm_insert_epi64(_mm_castpd_si128(force_x), _mm_cvtsd_f64(sp_x)*(fr * powTen), 0) );
            force_y = _mm_castsi128_pd (_mm_insert_epi64(_mm_castpd_si128(force_y), _mm_cvtsd_f64(sp_y)*(fr * powTen), 0) );
            force_z = _mm_castsi128_pd (_mm_insert_epi64(_mm_castpd_si128(force_z), _mm_cvtsd_f64(sp_z)*(fr * powTen), 0) );

        }
        
        if (rsqd1 >= 0.001 && rsqd1 < ptpCutOffSqd) {
            rsqd1 = rsqd1 * powTwenty;
            r = sqrt(rsqd1);
            rTwelve = rSix * rSix;
            rSix = ((double)rsqd1) * rsqd1 * rsqd1;
            f = (double)(VDW_A / rTwelve - VDW_B / rSix);
            if(doEnergy){
                energy += (double)( VDW_A / (12*rTwelve) - VDW_B / (6*rSix));
            }
            fr = f /r;
         
            force_x = _mm_castsi128_pd (_mm_insert_epi64(_mm_castpd_si128(force_x), _mm_cvtsd_f64(sp_x)*(fr * powTen), 1) );
            force_y = _mm_castsi128_pd (_mm_insert_epi64(_mm_castpd_si128(force_y), _mm_cvtsd_f64(sp_y)*(fr * powTen), 1) );
            force_z = _mm_castsi128_pd (_mm_insert_epi64(_mm_castpd_si128(force_z), _mm_cvtsd_f64(sp_z)*(fr * powTen), 1) );

         
        }
 
        fm_x = _mm_add_pd(fm_x, force_x);
        fm_y = _mm_add_pd(fm_y, force_y);
        fm_z = _mm_add_pd(fm_z, force_z);
        sm_x = _mm_add_pd(sm_x, force_x);
        sm_y = _mm_add_pd(sm_y, force_y);
        sm_z = _mm_add_pd(sm_z, force_z);
 
        //load the vectors back to the arrays
        _mm_store_pd(&first_x[i], f_x);
        _mm_store_pd(&first_y[i], f_y);
        _mm_store_pd(&first_z[i], f_z);
        _mm_store_pd(&sec_x[jpart], s_x);
        _mm_store_pd(&sec_y[jpart], s_y);
        _mm_store_pd(&sec_z[jpart], s_z);
        _mm_store_pd(&firstmsg_x[i], fm_x);
        _mm_store_pd(&firstmsg_y[i], fm_y);
        _mm_store_pd(&firstmsg_z[i], fm_z);
        _mm_store_pd(&secmsg_x[jpart], sm_x);
        _mm_store_pd(&secmsg_y[jpart], sm_y);
        _mm_store_pd(&secmsg_z[jpart], sm_z);
            
        
        }
    }
    
    
//copy back to original locations form the temp arrays
    
    //copy the values to the temp arrays
    for(int i=0; i<firstLen; i++){
        (first->part[i]).x = first_x[i];
        (first->part[i]).y = first_y[i];
        (first->part[i]).z = first_z[i];
        firstmsg[i].x = firstmsg_x[i];
        firstmsg[i].y = firstmsg_y[i];
        firstmsg[i].z = firstmsg_z[i];
        
    }
    for(int i=0; i<secondLen; i++){
        (second->part[i]).x = sec_x[i];
        (second->part[i]).y = sec_y[i];
        (second->part[i]).z = sec_z[i];
        secondmsg[i].x = secmsg_x[i];
        secondmsg[i].y = secmsg_y[i];
        secondmsg[i].z = secmsg_z[i];
    }
    
//sse end
    
/*
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
*/
    
  cellArray(first->x, first->y, first->z).receiveForces(stepCount,firstmsg,firstLen);
  cellArray(second->x, second->y, second->z).receiveForces(stepCount,secondmsg,secondLen);

  delete firstmsg;
  delete secondmsg;
  delete first;
  delete second;
  return energy;
}

//function to calculate forces among atoms in a single list
inline double calcInternalForces(ParticleDataMsg* first, int stepCount) {
  int i, j, ptpCutOffSqd;
  int firstLen = first->lengthAll;
  double powTwenty, powTen, firstx, firsty, firstz, rx, ry, rz, r, rsqd, fx, fy, fz, f, fr;
  vec3 firstpos, separation, force;
  double rSix, rTwelve;
  double energy = 0;
  int doEnergy = 0;
  if(stepCount == 1 || stepCount == finalStepCount)
    doEnergy = 1;
  vec3 *firstmsg = new vec3[firstLen];

  ptpCutOffSqd = PTP_CUT_OFF * PTP_CUT_OFF;
  powTen = pow(10.0, -10);
  powTwenty = pow(10.0, -20);

    //sse begin
    double  __attribute__ ((aligned (16))) first_x[firstLen], first_y[firstLen], first_z[firstLen];
    double  __attribute__ ((aligned (16))) firstmsg_x[firstLen], firstmsg_y[firstLen], firstmsg_z[firstLen];

    
    //copy the values to the temp arrays
    for(int i=0; i<firstLen; i++){
        first_x[i] = (first->part[i]).x;
        first_y[i] = (first->part[i]).y;
        first_z[i] = (first->part[i]).z;
        firstmsg_x[i] = firstmsg[i].x;
        firstmsg_y[i] = firstmsg[i].y;
        firstmsg_z[i] = firstmsg[i].z;
        
    }
    __m128d f_x, f_y, f_z, firstpos_x, firstpos_y, firstpos_z;
    __m128d fmi_x, fmi_y, fmi_z, fmj_x, fmj_y, fmj_z;
    __m128d sp_x, sp_y, sp_z, rsqd_x, rsqd_y, rsqd_z;
    __m128d force_x, force_y, force_z;
    
    for(i = 0; i < firstLen-3; i+=2){
        
        firstpos_x = _mm_load_pd(&first_x[i]);
		firstpos_y = _mm_load_pd(&first_y[i]);
		firstpos_z = _mm_load_pd(&first_z[i]);

        for(j = i; j < firstLen-3; j+=2) {
       //     ckout<<"x: "<< first_x[j]<< first_x[j+1]<< first_x[j+2]<< first_x[j+3]<<endl;
            f_x = _mm_load_pd(&first_x[j]);
            f_y = _mm_load_pd(&first_y[j]);
            f_z = _mm_load_pd(&first_z[j]);

            fmj_x = _mm_load_pd(&firstmsg_x[j]);
            fmj_y = _mm_load_pd(&firstmsg_y[j]);
            fmj_z = _mm_load_pd(&firstmsg_z[j]);
                                      
            fmi_x = _mm_load_pd(&firstmsg_x[i]);
            fmi_y = _mm_load_pd(&firstmsg_x[i]);
            fmi_z = _mm_load_pd(&firstmsg_x[i]);
            
            sp_x = _mm_sub_pd(firstpos_x, f_x); sp_y = _mm_sub_pd(firstpos_y, f_y); sp_z = _mm_sub_pd(firstpos_z, f_z);
            rsqd_x = _mm_mul_pd(sp_x, sp_x); rsqd_y = _mm_mul_pd(sp_y, sp_y); rsqd_z = _mm_mul_pd(sp_z, sp_z); //dot

            double rsqd0, rsqd1, rsqd2, rsqd3;
            rsqd0 = _mm_cvtsd_f64(rsqd_x) + _mm_cvtsd_f64(rsqd_y) + _mm_cvtsd_f64(rsqd_z);
            rsqd1 = _mm_cvtsd_f64(rsqd_x) + _mm_cvtsd_f64(rsqd_y) + _mm_cvtsd_f64(rsqd_z);

            force_x = _mm_set1_pd(0.0); force_y = _mm_set1_pd(0.0); force_z = _mm_set1_pd(0.0);
                                       

            if (rsqd0 >= 0.001 && rsqd0 < ptpCutOffSqd) {
                rsqd0 = rsqd0 * powTwenty;
                r = sqrt(rsqd0);
                rTwelve = rSix * rSix;
                rSix = ((double)rsqd0) * rsqd0 * rsqd0;
                f = (double)(VDW_A / rTwelve - VDW_B / rSix);
                if(doEnergy){
                    energy += (double)( VDW_A / (12*rTwelve) - VDW_B / (6*rSix));
                }
                fr = f /r;
                
                force_x = _mm_castsi128_pd (_mm_insert_epi32(_mm_castpd_si128(force_x), _mm_cvtsd_f64(sp_x)*(fr * powTen), 0) );
                force_y = _mm_castsi128_pd (_mm_insert_epi32(_mm_castpd_si128(force_y), _mm_cvtsd_f64(sp_y)*(fr * powTen), 0) );
                force_z = _mm_castsi128_pd (_mm_insert_epi32(_mm_castpd_si128(force_z), _mm_cvtsd_f64(sp_z)*(fr * powTen), 0) );
                
            }

            if (rsqd1 >= 0.001 && rsqd1 < ptpCutOffSqd) {
                rsqd1 = rsqd1 * powTwenty;
                r = sqrt(rsqd1);
                rTwelve = rSix * rSix;
                rSix = ((double)rsqd1) * rsqd1 * rsqd1;
                f = (double)(VDW_A / rTwelve - VDW_B / rSix);
                if(doEnergy){
                    energy += (double)( VDW_A / (12*rTwelve) - VDW_B / (6*rSix));
                }
                fr = f /r;
                force_x = _mm_castsi128_pd (_mm_insert_epi32(_mm_castpd_si128(force_x), _mm_cvtsd_f64(sp_x)*(fr * powTen), 1) );
                force_y = _mm_castsi128_pd (_mm_insert_epi32(_mm_castpd_si128(force_y), _mm_cvtsd_f64(sp_y)*(fr * powTen), 1) );
                force_z = _mm_castsi128_pd (_mm_insert_epi32(_mm_castpd_si128(force_z), _mm_cvtsd_f64(sp_z)*(fr * powTen), 1) );

                                          
                }
 
                fmj_x = _mm_add_pd(fmj_x, force_x);
                fmj_y = _mm_add_pd(fmj_y, force_y);
                fmj_z = _mm_add_pd(fmj_z, force_z);
                fmi_x = _mm_sub_pd(fmi_x, force_x);
                fmi_y = _mm_sub_pd(fmi_y, force_y);
                fmi_z = _mm_sub_pd(fmi_z, force_z);
 
                //load the vectors back to the arrays
                _mm_store_pd(&first_x[i], f_x);
                _mm_store_pd(&first_y[i], f_y);
                _mm_store_pd(&first_z[i], f_z);
                _mm_store_pd(&firstmsg_x[i], fmi_x);
                _mm_store_pd(&firstmsg_y[i], fmi_y);
                _mm_store_pd(&firstmsg_z[i], fmi_z);
                _mm_store_pd(&firstmsg_x[j], fmj_x);
                _mm_store_pd(&firstmsg_y[j], fmj_y);
                _mm_store_pd(&firstmsg_z[j], fmj_z);

        }
        _mm_store_pd(&first_x[i], firstpos_x);
        _mm_store_pd(&first_y[i], firstpos_y);
        _mm_store_pd(&first_z[i], firstpos_z);
                                      
    }
    
    //copy back to original locations form the temp arrays
                                      
    for(int i=0; i<firstLen; i++){
        (first->part[i]).x = first_x[i];
        (first->part[i]).y = first_y[i];
        (first->part[i]).z = first_z[i];
        firstmsg[i].x = firstmsg_x[i];
        firstmsg[i].y = firstmsg_y[i];
        firstmsg[i].z = firstmsg_z[i];
    }
    //sse end
   /*

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
    
    */
    
    
  cellArray(first->x, first->y, first->z).receiveForces(stepCount,firstmsg,firstLen);
  delete firstmsg;
  delete first;
  return energy;
}
#endif
