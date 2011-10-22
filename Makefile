#for compilation do a make

#for run use following - Either provide no arguments or provide all
#./charmrun +p<Procs> ./leanmd <dimX,dimY,dimZ,steps,firstLBstep,LBPeriod,ftPeriod>
# Example - ./charmrun +p8 ./leanmd 2 2 2 1001 20 20 10000

CHARMBASE	= $(HOME)/charm/bluegenep-xlc
CHARMC          = $(CHARMBASE)/bin/charmc
OPTS            = -O3 

all: leanmd

leanmd: Main.o Patch.o Compute.o leanmd.decl.h
	$(CHARMC) $(OPTS) -module CkMulticast -module CommonLBs \
	-module HybridLB -language charm++ -o leanmd Main.o Patch.o Compute.o

Main.o: Main.C Main.h leanmd.decl.h defs.h
	$(CHARMC) $(OPTS) -o Main.o Main.C

Patch.o: Patch.C Patch.h leanmd.decl.h defs.h
	$(CHARMC) $(OPTS) -o Patch.o Patch.C

leanmd.decl.h:	leanmd.ci
	$(CHARMC) leanmd.ci

Compute.o: Compute.C Compute.h leanmd.decl.h defs.h physics.h
	$(CHARMC) $(OPTS) -o Compute.o Compute.C


clean:
	rm -f *.decl.h *.def.h *.o leanmd leanmd.prj charmrun
