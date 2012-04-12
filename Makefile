# to be set appropiately
CHARMBASE      = /expand/home/nikhil/charms/charm/net-linux-x86_64
CHARMC         ?= $(CHARMBASE)/bin/charmc

OPTS            = -O3 -g

all: leanmd

leanmd: Main.o Cell.o Compute.o leanmd.decl.h
	$(CHARMC) $(OPTS) -tracemode projections -module CkMulticast -module CommonLBs -language charm++ -o leanmd Main.o Cell.o Compute.o

Main.o: Main.cc Main.h Comm.h leanmd.decl.h defs.h
	$(CHARMC) $(OPTS) -o Main.o Main.cc

Cell.o: Cell.cc Cell.h Comm.h leanmd.decl.h defs.h
	$(CHARMC) $(OPTS) -o Cell.o Cell.cc

leanmd.decl.h:	leanmd.ci Comm.h
	$(CHARMC) leanmd.ci

Compute.o: Compute.cc Compute.h Comm.h leanmd.decl.h defs.h physics.h
	$(CHARMC) $(OPTS) -o Compute.o Compute.cc

test: leanmd
	./charmrun +p4 ./leanmd 4 4 4 101 20 20 +balancer GreedyLB +LBDebug 1

clean:
	rm -f *.decl.h *.def.h *.o leanmd leanmd.prj charmrun
