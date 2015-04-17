SHELL := /bin/bash

# to be set appropiately
CHARMBASE	=	$(HOME)/PPL/git/charmgitrrit/charm/mpi-linux-x86_64
#CHARMBASE      = $(HOME)/PPL/git/charmSteering/charm/mpi-linux-x86_64
#CHARMBASE      = $(HOME)/git/charmSteering/charm/pamilrts-bluegeneq
CHARMC         = $(CHARMBASE)/bin/charmc

OPTS            = -O3 -L$(HOME)/papi/lib -g

DECL=
SUFFIX=

ifeq ($(FT), yes)
DECL=-DCMK_MEM_CHECKPOINT=1
FT=-syncft
endif

all: leanmd

leanmd: Main.o Cell.o Compute.o leanmd.decl.h
	$(CHARMC) $(OPTS) -module CkMulticast -module CommonLBs \
	-language charm++ -o leanmd$(SUFFIX) Main.o Cell.o Compute.o -tracemode autoPerf -tracemode autoTuner

leanmd.prj: Main.o Cell.o Compute.o leanmd.decl.h
	$(CHARMC) $(OPTS) -module CkMulticast -module CommonLBs \
	-language charm++ -o leanmd.prj Main.o Cell.o Compute.o -tracemode projections  

Main.o: Main.cc Main.h leanmd.decl.h defs.h
	$(CHARMC) $(OPTS) -o Main.o Main.cc

Cell.o: Cell.cc Cell.h leanmd.decl.h defs.h
	$(CHARMC) $(OPTS) -o Cell.o Cell.cc

leanmd.decl.h:	leanmd.ci
	$(CHARMC) -E leanmd.ci $(DECL)

Compute.o: Compute.cc Compute.h leanmd.decl.h defs.h physics.h
	$(CHARMC) $(OPTS) -o Compute.o Compute.cc

test: leanmd
	./charmrun +p4 ./leanmd 4 4 4 10 3 3 10000 mem 1 1 1  4 +balancer GreedyLB +LBDebug 1 ++local

clean:
	rm -f *.decl.h *.def.h *.o leanmd leanmd-ft leanmd.prj charmrun
