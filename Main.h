#ifndef __MAIN_H__
#define __MAIN_H__

//central controller chare
class Main : public CBase_Main {
  private:
    int phase;    //variable to keep track of phase in initial set up
    double energy, prevEnergy;  //store initial and final energy
    int testFailed;    //flag to indicate if simulation conserve energy
    int endCount;

  public:
    Main(CkArgMsg* msg);
    Main(CkMigrateMessage* msg);
    void pup(PUP::er &p);

    //member functions
    void allDone();
    void ftBarrier();
    void startUpDone();
    void energySum(double energyIn);
};
#endif
