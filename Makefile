# to be set appropiately
CHARMBASE      = /expand/home/nikhil/charms/charm/net-linux-x86_64
CHARMC         = $(CHARMBASE)/bin/charmc

CXXFLAGS            = -O0 -g
CXX = $(CHARMC)

all: leanmd

leanmd: Main.o Cell.o Compute.o ComputePME.o
	$(CHARMC) $(OPTS) -module CkMulticast -module CommonLBs \
	-language charm++ -o leanmd $^

Main.o: Main.cc Main.h leanmd.decl.h defs.h
Cell.o: Cell.cc Cell.h leanmd.decl.h defs.h
Compute.o: Compute.cc Compute.h leanmd.decl.h defs.h physics.h
ComputePME.o: ComputePME.cc ComputePME.h leanmd.decl.h defs.h Cell.h

leanmd.decl.h:	leanmd.ci
	$(CHARMC) leanmd.ci


test: leanmd
	./charmrun +p4 ./leanmd 4 4 4 101 20 20 +balancer GreedyLB +LBDebug 1

clean:
	rm -f *.decl.h *.def.h *.o leanmd leanmd.prj charmrun
