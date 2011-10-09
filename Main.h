#ifndef __MAIN_H__
#define __MAIN_H__

class Main : public CBase_Main {
  private:
    int phase;

  public:
    Main(CkArgMsg* msg);
    Main(CkMigrateMessage* msg);
    void pup(PUP::er &p) {
      Chare::pup(p);
      p|phase;
    }

    void allDone();
    void lbBarrier();
    void ftBarrier();
    void startUpDone();
};
#endif
