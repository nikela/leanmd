#ifndef __MAIN_H__
#define __MAIN_H__

extern /* readonly */ CProxy_Patch patchArray;
extern /* readonly */ CProxy_Compute computeArray;

extern /* readonly */ int patchArrayDimX;
extern /* readonly */ int patchArrayDimY;
extern /* readonly */ int patchArrayDimZ;

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

    //member functions
    //void allDone();
    //void energySum(double energyIn);
};
#endif
