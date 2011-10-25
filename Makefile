include ../config.mk

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

test:
	./charmrun +p4 ./leanmd 4 4 4 101 20 20

clean:
	rm -f *.decl.h *.def.h *.o leanmd leanmd.prj charmrun
