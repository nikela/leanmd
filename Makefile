SHELL := /bin/bash

# to be set appropiately
CHARMBASE      = $(HOME)/charm/
CHARMC         = $(CHARMBASE)/bin/charmc
ifeq ($(TARGET),gni)
  MPICXX=CC
else
ifeq ($(TARGET),bgq)
  MPICXX=mpicxx
  BGQ_INSTALL=/bgsys/drivers/ppcfloor
  MPI_LIBS=-L$(BGQ_INSTALL)/comm/gcc/lib -lmpich -lopa -lmpl -ldl -L$(BGQ_INSTALL)/comm/sys-fast/lib -lpami -L$(BGQ_INSTALL)/spi/lib -lSPI -lSPI_cnk -lpthread -lrt
else
ifeq ($(TARGET),xlc)
  MPICXX=mpixlcxx_r
  BGQ_INSTALL=/bgsys/drivers/ppcfloor
  MPI_LIBS=-L$(BGQ_INSTALL)/comm/xl/lib -lmpich -lopa -lmpl -ldl -L$(BGQ_INSTALL)/comm/sys-fast/lib -lpami -L$(BGQ_INSTALL)/spi/lib -lSPI -lSPI_cnk -lpthread -lrt
else
  MPICXX=mpicxx
endif
endif
endif

OPTS            = -O0 -g

DECL=
SUFFIX=

ifeq ($(FT), yes)
DECL=-DCMK_MEM_CHECKPOINT=1
FT=-syncft
endif

all: interop

interop: leanmd LeanMD.cc
	$(MPICXX) $(OPTS) -c LeanMD.cc -I$(CHARMBASE)/include
	$(CHARMC) $(OPTS) -mpi -nomain-module -language charm++ -o leanmd-interop LeanMD.o -L. -module leanmd -module CkMulticast -module CommonLBs -module ParMetisLB -module ParMetisCentLB -lparmetis -lmetis $(MPI_LIBS)

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
