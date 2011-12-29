#ifndef __MAIN_H__
#define __MAIN_H__

extern /* readonly */ CProxy_Cell cellArray;
extern /* readonly */ CProxy_Compute computeArray;

extern /* readonly */ int cellArrayDimX;
extern /* readonly */ int cellArrayDimY;
extern /* readonly */ int cellArrayDimZ;

//central controller chare
class Main : public CBase_Main {
  private:
    Main_SDAG_CODE
    double energy, prevEnergy;  //store initial and final energy
    int testFailed;    //flag to indicate if simulation conserve energy
    int endCount;

  public:
    Main(CkArgMsg* msg);
    Main(CkMigrateMessage* msg);
    void pup(PUP::er &p);
};
#endif
