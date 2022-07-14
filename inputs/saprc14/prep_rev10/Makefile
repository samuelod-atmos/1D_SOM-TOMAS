# Generated automatically from Makefile.in by configure.
ALL: default
##### User configurable options #####

SHELL       = /bin/sh
F77         = ifort
FLINKER     = $(F77)
#

### End User configurable options ###

#FFLAGS	  = $(OPTFLAGS)
#FFLAGS    = -Mbounds
FFLAGS = -static
#FFLAGS    = -Vaxlib
EXECS	  = prep.exe prep_bag.exe

default: $(EXECS)

all: default

prep.exe: prep.o asmd.o bldup.o diffun.o inscod.o listsc.o newinbuf.o \
          newsubs.o pfalin.o pnhrxn.o rdrxn.o typeo.o reordr.o rxlst1.o \
          rxlst2.o rxsubs.o saverx.o typeosubs.o vcoef.o som.o somsubs.o
	$(F77) $(FFLAGS) -o prep.exe \
       prep.o asmd.o bldup.o diffun.o inscod.o listsc.o newinbuf.o \
       newsubs.o pfalin.o pnhrxn.o rdrxn.o typeo.o reordr.o rxlst1.o \
       rxlst2.o rxsubs.o saverx.o typeosubs.o vcoef.o som.o somsubs.o

prep_bag.exe: prep.o asmd.o bldup.o diffun.o inscod.o listsc.o newinbuf.o \
          newsubs.o pfalin.o pnhrxn.o rdrxn.o typeo.o reordr.o rxlst1.o \
          rxlst2.o rxsubs.o saverx_bag.o typeosubs.o vcoef.o som.o somsubs.o
	$(F77) $(FFLAGS) -o prep_bag.exe \
       prep.o asmd.o bldup.o diffun.o inscod.o listsc.o newinbuf.o \
       newsubs.o pfalin.o pnhrxn.o rdrxn.o typeo.o reordr.o rxlst1.o \
       rxlst2.o rxsubs.o saverx_bag.o typeosubs.o vcoef.o som.o somsubs.o

prep.o: prep.for *.inc
	$(F77) $(FFLAGS) -c prep.for

asmd.o: asmd.for *.inc
	$(F77) $(FFLAGS) -c asmd.for

bldup.o: bldup.for *.inc
	$(F77) $(FFLAGS) -c bldup.for

diffun.o: diffun.for *.inc
	$(F77) $(FFLAGS) -c diffun.for

inscod.o: inscod.for *.inc
	$(F77) $(FFLAGS) -c inscod.for

listsc.o: listsc.for *.inc
	$(F77) $(FFLAGS) -c listsc.for

newinbuf.o: newinbuf.for *.inc
	$(F77) $(FFLAGS) -c newinbuf.for

newsubs.o: newsubs.for *.inc
	$(F77) $(FFLAGS) -c newsubs.for 

pfalin.o: pfalin.for *.inc
	$(F77) $(FFLAGS) -c pfalin.for

pnhrxn.o: pnhrxn.for *.inc
	$(F77) $(FFLAGS) -c pnhrxn.for

rdrxn.o: rdrxn.for *.inc
	$(F77) $(FFLAGS) -c rdrxn.for

typeo.o: typeo.for *.inc
	$(F77) $(FFLAGS) -c typeo.for

typeosubs.o: typeosubs.for *.inc
	$(F77) $(FFLAGS) -c typeosubs.for

reordr.o: reordr.for *.inc
	$(F77) $(FFLAGS) -c reordr.for

rxlst1.o: rxlst1.for *.inc
	$(F77) $(FFLAGS) -c rxlst1.for

rxlst2.o: rxlst2.for *.inc
	$(F77) $(FFLAGS) -c rxlst2.for

rxsubs.o: rxsubs.for *.inc
	$(F77) $(FFLAGS) -c rxsubs.for

saverx.o: saverx.for *.inc
	$(F77) $(FFLAGS) -c saverx.for

saverx_bag.o: saverx_bag.for *.inc
	$(F77) $(FFLAGS) -c saverx_bag.for

vcoef.o: vcoef.for *.inc
	$(F77) $(FFLAGS) -c vcoef.for

som.o: som.for *.inc
	$(F77) $(FFLAGS) -c som.for

somsubs.o: somsubs.for *.inc
	$(F77) $(FFLAGS) -c somsubs.for

clean:
	/bin/rm -f *.o *~ PI* $(EXECS) rdb.* startup.*

