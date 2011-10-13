CHARMBASE	= $(HOME)/codes/charm/net-linux-x86_64
CHARMC          = $(CHARMBASE)/bin/charmc
OPTS            = -O3 -DUSE_SECTION_MULTICAST

all: leanmd

leanmd: Main.o Patch.o Compute.o leanmd.decl.h
	$(CHARMC) $(OPTS) -module CkMulticast -module CommonLBs \
	-language charm++ -o leanmd Main.o Patch.o Compute.o

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
