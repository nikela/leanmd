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
    double finalEnergy, initialEnergy;  //store initial and final energy
    int testFailed;    //flag to indicate if the simulation conserves energy

  public:
    Main(CkArgMsg* msg);
    Main(CkMigrateMessage* msg);
    void pup(PUP::er &p);
};
#endif
