
#ifndef __DEFS__
#define __DEFS__

#include "pup.h"

#define HYDROGEN_MASS   (1.67 * pow( 10.0,-24))
#define VDW_A			      (1.60694452 * pow(10.0, -134))
#define VDW_B			      (1.031093844 * pow(10.0, -77))

#define ENERGY_VAR  		(1.0 * pow(10.0,-5))
#define PARTICLES_PER_PATCH   	150	

#define DEFAULT_DELTA		1	// in femtoseconds

#define DEFAULT_FIRST_LDB	20
#define DEFAULT_LDB_PERIOD	20
#define DEFAULT_FT_PERIOD	100000

#define PATCHARRAY_DIM_X	3
#define PATCHARRAY_DIM_Y	3
#define PATCHARRAY_DIM_Z	3
#define PTP_CUT_OFF		12	//  cut off for atom to atom interactions
#define PATCH_MARGIN		4 	// constant difference between cut off and patch size
#define PATCH_SIZE_X		(PTP_CUT_OFF + PATCH_MARGIN)
#define PATCH_SIZE_Y		(PTP_CUT_OFF + PATCH_MARGIN)
#define PATCH_SIZE_Z		(PTP_CUT_OFF + PATCH_MARGIN)
#define PATCH_ORIGIN_X		0
#define PATCH_ORIGIN_Y		0
#define PATCH_ORIGIN_Z		0

#define MIGRATE_STEPCOUNT	20
#define DEFAULT_FINALSTEPCOUNT	1001
#define MAX_VELOCITY		30.0

#define KAWAY_X			2
#define KAWAY_Y			2
#define KAWAY_Z			1
#define NBRS_X			(2*KAWAY_X+1)
#define NBRS_Y			(2*KAWAY_Y+1)
#define NBRS_Z			(2*KAWAY_Z+1)
#define NUM_NEIGHBORS		(NBRS_X * NBRS_Y * NBRS_Z)

#define WRAP_X(a)		(((a)+patchArrayDimX)%patchArrayDimX)
#define WRAP_Y(a)		(((a)+patchArrayDimY)%patchArrayDimY)
#define WRAP_Z(a)		(((a)+patchArrayDimZ)%patchArrayDimZ)

struct vec3 {
  double x, y, z;

  vec3(double d = 0.0) : x(d), y(d), z(d) { }
  vec3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) { }

  inline vec3& operator += (const vec3 &rhs) {
    x += rhs.x; y += rhs.y; z += rhs.z;
    return *this;
  }
  inline vec3& operator -= (const vec3 &rhs) {
    return *this += (rhs * -1.0);
  }
  inline vec3 operator* (const double d) const {
    return vec3(d*x, d*y, d*z);
  }
  inline vec3 operator- (const vec3& rhs) const {
    return vec3(x - rhs.x, y - rhs.y, z - rhs.z);
  }
};
inline double dot(const vec3& a, const vec3& b) {
  return a.x*b.x + a.y*b.y + a.z*b.z;
}
PUPbytes(vec3)

//class for keeping track of the properties for a particle
struct Particle {
//  int id;
  double mass;
  //   Position, acceleration, velocity
  vec3 pos,      acc,          vel;

  // Function for pupping properties
  void pup(PUP::er &p) {
    p | mass;
    p | pos;
    p | acc;
    p | vel;
  }
};
#endif
