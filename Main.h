#ifndef __MAIN_H__
#define __MAIN_H__

//central controller chare
class Main : public CBase_Main {
  private:
    int phase;    //variable to keep track of phase in initial set up
    double energy, prevEnergy;  //store initial and final energy
    int testFailed;    //flag to indicate if simulation conserve energy

  public:
    Main(CkArgMsg* msg);
    Main(CkMigrateMessage* msg);

    //pup routine incase the main chare moves, pack important information
    void pup(PUP::er &p) {
      Chare::pup(p);
      p|phase;
      p|energy;
      p|prevEnergy;
      p|testFailed;
    }

    //member functions
    void allDone();
    void ftBarrier();
    void startUpDone();
    void energySumP(double energyP);
    void energySumK(double energyK);
};
#endif
