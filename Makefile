#for compilation do a make

#for run use following - Either provide no arguments or provide all
#./charmrun +p<Procs> ./leanmd <dimX,dimY,dimZ,steps,firstLBstep,LBPeriod>
# Example - ./charmrun +p8 ./leanmd 2 2 2 1001 20 20 

CHARMBASE				=	$(HOME)/codes/charm/net-linux-x86_64
CHARMC          = $(CHARMBASE)/bin/charmc
OPTS            = -O3 

all: leanmd

leanmd: Main.o Patch.o Compute.o leanmd.decl.h
	$(CHARMC) $(OPTS) -module CkMulticast -module CommonLBs \
	-module HybridLB -language charm++ -o leanmd Main.o Patch.o Compute.o

Main.o: Main.cc Main.h leanmd.decl.h defs.h
	$(CHARMC) $(OPTS) -o Main.o Main.cc

Patch.o: Patch.cc Patch.h leanmd.decl.h defs.h
	$(CHARMC) $(OPTS) -o Patch.o Patch.cc

leanmd.decl.h:	leanmd.ci
	$(CHARMC) leanmd.ci

Compute.o: Compute.cc Compute.h leanmd.decl.h defs.h physics.h
	$(CHARMC) $(OPTS) -o Compute.o Compute.cc


clean:
	rm -f *.decl.h *.def.h *.o leanmd leanmd.prj charmrun
