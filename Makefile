SHELL := /bin/bash

# to be set appropiately
CHARMBASE      = $(HOME)/PPL/git/charmSteering/charm/mpi-linux-x86_64
CHARMC         = $(CHARMBASE)/bin/charmc

OPTS            = -O3 -L/dcsdata/home/jessie/papi/lib

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
	-language charm++ -o leanmd.prj Main.o Cell.o Compute.o -tracemode projections  -tracemode autoPerf -tracemode autoTuner

Main.o: Main.cc Main.h leanmd.decl.h defs.h
	$(CHARMC) $(OPTS) -o Main.o Main.cc

Cell.o: Cell.cc Cell.h leanmd.decl.h defs.h
	$(CHARMC) $(OPTS) -o Cell.o Cell.cc

leanmd.decl.h:	leanmd.ci
	$(CHARMC) -E leanmd.ci $(DECL)

Compute.o: Compute.cc Compute.h leanmd.decl.h defs.h physics.h
	$(CHARMC) $(OPTS) -o Compute.o Compute.cc

test: leanmd
	./charmrun +p4 ./leanmd 4 4 4 10 3 3 +balancer GreedyLB +LBDebug 1 ++local

clean:
	rm -f *.decl.h *.def.h *.o leanmd leanmd-ft leanmd.prj charmrun
