#ifndef __PHYSICS_H__
#define __PHYSICS_H__


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

//qpx added

double __attribute__((aligned(32))) first_x[firstLen], first_y[firstLen], first_z[firstLen];
double __attribute__((aligned(32))) sec_x[firstLen], sec_y[firstLen], sec_z[firstLen];
double __attribute__((aligned(32))) firstmsg_x[firstLen], firstmsg_y[firstLen], firstmsg_z[firstLen];
double __attribute__((aligned(32))) secmsg_x[firstLen], secmsg_y[firstLen], secmsg_z[firstLen];

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
    
vector4double f_x, f_y, f_z, s_x, s_y, s_z;
vector4double fm_x, fm_y, fm_z, sm_x, sm_y, sm_z;
vector4double sp_x, sp_y, sp_z, rsqd_x, rsqd_y, rsqd_z;
vector4double force_x, force_y, force_z;
    
  int i1, j1;
  for(i1 = 0; i1 < firstLen; i1=i1+BLOCK_SIZE)
    for(j1 = 0; j1 < secondLen; j1=j1+BLOCK_SIZE)
      for(i = i1; i < i1+BLOCK_SIZE && i < firstLen; i+=4) {
        for(jpart = j1; jpart < j1+BLOCK_SIZE && jpart < secondLen; jpart+=4) {

        //load the arrays to vectors
		f_x = vec_ld(0L, &first_x[i]);
		f_y = vec_ld(0L, &first_y[i]);
		f_z = vec_ld(0L, &first_z[i]);		
		s_x = vec_ld(0L, &sec_x[jpart]);
		s_y = vec_ld(0L, &sec_y[jpart]);
		s_z = vec_ld(0L, &sec_z[jpart]);
        fm_x = vec_ld(0L, &firstmsg_x[i]);
        fm_y = vec_ld(0L, &firstmsg_y[i]);
        fm_z = vec_ld(0L, &firstmsg_z[i]);
        sm_x = vec_ld(0L, &secmsg_x[jpart]);
        sm_y = vec_ld(0L, &secmsg_y[jpart]);
        sm_z = vec_ld(0L, &secmsg_z[jpart]);
            
        sp_x = vec_sub(f_x, s_x); sp_y = vec_sub(f_y, s_y); sp_z = vec_sub(f_z, s_z);
        rsqd_x = vec_mul(sp_x, sp_x); rsqd_y = vec_mul(sp_y, sp_y); rsqd_z = vec_mul(sp_z, sp_z); //dot
        
        double rsqd0, rsqd1, rsqd2, rsqd3;
        rqsd0 = vec_extract(rsqd_x, 0) + vec_extract(rsqd_y, 0) + vec_extract(rsqd_z, 0);
        rqsd1 = vec_extract(rsqd_x, 1) + vec_extract(rsqd_y, 1) + vec_extract(rsqd_z, 1);
        rqsd2 = vec_extract(rsqd_x, 2) + vec_extract(rsqd_y, 2) + vec_extract(rsqd_z, 2);
		rqsd3 = vec_extract(rsqd_x, 3) + vec_extract(rsqd_y, 3) + vec_extract(rsqd_z, 3);
        force_x = vec_splats(0.0); force_y = vec_splats(0.0); force_z = vec_splats(0.0);
        
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
            force_x = vec_insert( vec_extract(sp_x, 0) * (fr * powTen), force_x, 0 );
            force_y = vec_insert( vec_extract(sp_y, 0) * (fr * powTen), force_y, 0 );
            force_z = vec_insert( vec_extract(sp_z, 0) * (fr * powTen), force_z, 0 );
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
            force_x = vec_insert( vec_extract(sp_x, 1) * (fr * powTen), force_x, 1 );
            force_y = vec_insert( vec_extract(sp_y, 1) * (fr * powTen), force_y, 1 );
            force_z = vec_insert( vec_extract(sp_z, 1) * (fr * powTen), force_z, 1 );
                
        }
        if (rsqd2 >= 0.001 && rsqd2 < ptpCutOffSqd) {
            rsqd2 = rsqd2 * powTwenty;
            r = sqrt(rsqd2);
            rTwelve = rSix * rSix;
            rSix = ((double)rsqd2) * rsqd2 * rsqd2;
            f = (double)(VDW_A / rTwelve - VDW_B / rSix);
            if(doEnergy){
                energy += (double)( VDW_A / (12*rTwelve) - VDW_B / (6*rSix));
            }
            fr = f /r;
            force_x = vec_insert( vec_extract(sp_x, 2) * (fr * powTen), force_x, 2 );
            force_y = vec_insert( vec_extract(sp_y, 2) * (fr * powTen), force_y, 2 );
            force_z = vec_insert( vec_extract(sp_z, 2) * (fr * powTen), force_z, 2 );
                
        }
        if (rsqd3 >= 0.001 && rsqd3 < ptpCutOffSqd) {
            rsqd3 = rsqd3 * powTwenty;
            rSix = ((double)rsqd3) * rsqd3 * rsqd3;
            r = sqrt(rsqd3);
            rTwelve = rSix * rSix;
            f = (double)(VDW_A / rTwelve - VDW_B / rSix);
            if(doEnergy){
                energy += (double)( VDW_A / (12*rTwelve) - VDW_B / (6*rSix));
            }
            fr = f /r;
            force_x = vec_insert( vec_extract(sp_x, 3) * (fr * powTen), force_x, 3 );
            force_y = vec_insert( vec_extract(sp_y, 3) * (fr * powTen), force_y, 3 );
            force_z = vec_insert( vec_extract(sp_z, 3) * (fr * powTen), force_z, 3 );
        }
        //add force to... sub force from ..
        fm_x = vec_add(fm_x, force_x);
        fm_y = vec_add(fm_y, force_y);
        fm_z = vec_add(fm_z, force_z);
        sm_x = vec_sub(sm_x, force_x);
        sm_y = vec_sub(sm_y, force_y);
        sm_z = vec_sub(sm_z, force_z);
        
        //load the vectors back to the arrays
        vec_st(f_x, 0L, &first_x[i]);
        vec_st(f_y, 0L, &first_y[i]);
        vec_st(f_z, 0L, &first_z[i]);
        vec_st(s_x, 0L, &sec_x[jpart]);
        vec_st(s_y, 0L, &sec_y[jpart]);
        vec_st(s_z, 0L, &sec_z[jpart]);
        vec_st(fm_x, 0L, &firstmsg_x[i]);
        vec_st(fm_y, 0L, &firstmsg_y[i]);
        vec_st(fm_z, 0L, &firstmsg_z[i]);
        vec_st(sm_x, 0L, &secmsg_x[jpart]);
        vec_st(sm_y, 0L, &secmsg_y[jpart]);
        vec_st(sm_z, 0L, &secmsg_z[jpart]);           
            
        
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

//qpx end

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

    //qpx added
 
    double __attribute__((aligned(32))) first_x[firstLen], first_y[firstLen], first_z[firstLen];
    double __attribute__((aligned(32))) firstmsg_x[firstLen], firstmsg_y[firstLen], firstmsg_z[firstLen];

    
    //copy the values to the temp arrays
    for(int i=0; i<firstLen; i++){
        first_x[i] = (first->part[i]).x;
        first_y[i] = (first->part[i]).y;
        first_z[i] = (first->part[i]).z;
        firstmsg_x[i] = firstmsg[i].x;
        firstmsg_y[i] = firstmsg[i].y;
        firstmsg_z[i] = firstmsg[i].z;
        
    }
    vector4double f_x, f_y, f_z, firstpos_x, firstpos_y, firstpos_z;
    vector4double fmi_x, fmi_y, fmi_z, fmj_x, fmj_y, fmj_z;
    vector4double sp_x, sp_y, sp_z, rsqd_x, rsqd_y, rsqd_z;
    vector4double force_x, force_y, force_z;
    
    for(i = 0; i < firstLen; i+=4){
        firstpos_x = vec_ld(0L, &first_x[i]);
		firstpos_y = vec_ld(0L, &first_y[i]);
		firstpos_z = vec_ld(0L, &first_z[i]);

        for(j = i+1; j < firstLen; j+=4) {
            
            f_x = vec_ld(0L, &first_x[j);
            f_y = vec_ld(0L, &first_y[j]);
            f_z = vec_ld(0L, &first_z[j]);
                                      
            fmj_x = vec_ld(0L, &firstmsg_x[j]);
            fmj_y = vec_ld(0L, &firstmsg_y[j]);
            fmj_z = vec_ld(0L, &firstmsg_z[j]);
                                      
            fmi_x = vec_ld(0L, &firstmsg_x[i]);
            fmi_y = vec_ld(0L, &firstmsg_x[i]);
            fmi_z = vec_ld(0L, &firstmsg_x[i]);
            
            sp_x = vec_sub(firstpos_x, f_x); sp_y = vec_sub(firstpos_y, f_y); sp_z = vec_sub(firstpos_z, f_z);
            rsqd_x = vec_mul(sp_x, sp_x); rsqd_y = vec_mul(sp_y, sp_y); rsqd_z = vec_mul(sp_z, sp_z); //dot
            
            double rsqd0, rsqd1, rsqd2, rsqd3;
            rqsd0 = vec_extract(rsqd_x, 0) + vec_extract(rsqd_y, 0) + vec_extract(rsqd_z, 0);
            rqsd1 = vec_extract(rsqd_x, 1) + vec_extract(rsqd_y, 1) + vec_extract(rsqd_z, 1);
            rqsd2 = vec_extract(rsqd_x, 2) + vec_extract(rsqd_y, 2) + vec_extract(rsqd_z, 2);
            rqsd3 = vec_extract(rsqd_x, 3) + vec_extract(rsqd_y, 3) + vec_extract(rsqd_z, 3);
            force_x = vec_splats(0.0); force_y = vec_splats(0.0); force_z = vec_splats(0.0);
                                      
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
                force_x = vec_insert( vec_extract(sp_x, 0) * (fr * powTen), force_x, 0 );
                force_y = vec_insert( vec_extract(sp_y, 0) * (fr * powTen), force_y, 0 );
                force_z = vec_insert( vec_extract(sp_z, 0) * (fr * powTen), force_z, 0 );
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
                force_x = vec_insert( vec_extract(sp_x, 1) * (fr * powTen), force_x, 1 );
                force_y = vec_insert( vec_extract(sp_y, 1) * (fr * powTen), force_y, 1 );
                force_z = vec_insert( vec_extract(sp_z, 1) * (fr * powTen), force_z, 1 );
                                          
                }
                if (rsqd2 >= 0.001 && rsqd2 < ptpCutOffSqd) {
                    rsqd2 = rsqd2 * powTwenty;
                    r = sqrt(rsqd2);
                    rTwelve = rSix * rSix;
                    rSix = ((double)rsqd2) * rsqd2 * rsqd2;
                    f = (double)(VDW_A / rTwelve - VDW_B / rSix);
                    if(doEnergy){
                        energy += (double)( VDW_A / (12*rTwelve) - VDW_B / (6*rSix));
                    }
                    fr = f /r;
                    force_x = vec_insert( vec_extract(sp_x, 2) * (fr * powTen), force_x, 2 );
                    force_y = vec_insert( vec_extract(sp_y, 2) * (fr * powTen), force_y, 2 );
                    force_z = vec_insert( vec_extract(sp_z, 2) * (fr * powTen), force_z, 2 );
                                          
                }
                if (rsqd3 >= 0.001 && rsqd3 < ptpCutOffSqd) {
                    rsqd3 = rsqd3 * powTwenty;
                    rSix = ((double)rsqd3) * rsqd3 * rsqd3;
                    r = sqrt(rsqd3);
                    rTwelve = rSix * rSix;
                    f = (double)(VDW_A / rTwelve - VDW_B / rSix);
                    if(doEnergy){
                        energy += (double)( VDW_A / (12*rTwelve) - VDW_B / (6*rSix));
                    }
                    fr = f /r;
                    force_x = vec_insert( vec_extract(sp_x, 3) * (fr * powTen), force_x, 3 );
                    force_y = vec_insert( vec_extract(sp_y, 3) * (fr * powTen), force_y, 3 );
                    force_z = vec_insert( vec_extract(sp_z, 3) * (fr * powTen), force_z, 3 );
                }
                //add force to... sub force from ..
                fmj_x = vec_add(fmj_x, force_x);
                fmj_y = vec_add(fmj_y, force_y);
                fmj_z = vec_add(fmj_z, force_z);
                fmi_x = vec_sub(fmi_x, force_x);
                fmi_y = vec_sub(fmi_y, force_y);
                fmi_z = vec_sub(fmi_z, force_z);
                                      
                //load the vectors back to the arrays
                vec_st(f_x, 0L, &first_x[j]);
                vec_st(f_y, 0L, &first_y[j]);
                vec_st(f_z, 0L, &first_z[j]);
                vec_st(fmi_x, 0L, &firstmsg_x[i]);
                vec_st(fmi_y, 0L, &firstmsg_y[i]);
                vec_st(fmi_z, 0L, &firstmsg_z[i]);
                vec_st(fmj_x, 0L, &firstmsg_x[j]);
                vec_st(fmj_y, 0L, &firstmsg_y[j]);
                vec_st(fmj_z, 0L, &firstmsg_z[j]);
        }
        vec_st(firstpos_x, 0L, &first_x[i]);
        vec_st(firstpos_y, 0L, &first_y[i]);
        vec_st(firstpos_z, 0L, &first_z[i]);
                                      
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
                                
                                      
    //qpx end
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
