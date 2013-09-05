SHELL := /bin/bash

# to be set appropiately
CHARMBASE      = $(HOME)/work/charm
CHARMC         = $(CHARMBASE)/bin/charmc
MPICXX         = mpicxx

OPTS            = -O3

DECL=
SUFFIX=

ifeq ($(FT), yes)
DECL=-DCMK_MEM_CHECKPOINT=1
FT=-syncft
endif

all: interop

interop: leanmd LeanMD.cc
	$(MPICXX) $(OPTS) -c LeanMD.cc -I$(CHARMBASE)/include
	$(CHARMC) $(OPTS) -mpi -nomain-module -language charm++ -o leanmd-interop LeanMD.o -L. -module leanmd -module CkMulticast -module CommonLBs -module HybridLB

leanmd: Main.o Cell.o Compute.o leanmd.decl.h
	$(CHARMC) $(OPTS) -language charm++ -o libmoduleleanmd.a Main.o Cell.o Compute.o

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
	rm -f *.decl.h *.def.h *.o leanmd leanmd-ft leanmd.prj libmoduleleanmd.a charmrun
