
# Makefile for TOMAS box model
# questions or comments, email Jeff Pierce, jrpierce@andrew.cmu.edu
###########################################################################

# options
comp = gfortran
# to run quickly, use -O2 flag (optimize level 2) or -O if using g77 with no -g flag (for debugger)
# to debug, remove -O2 and add -g flag
o_flags = -c -O1 -fcheck=all -fbounds-check #-v # flags for the compiling of *.o files [O = optimization] 
exe_flags = -O1 -dynamic -fcheck=all -fbounds-check # flags for the compiling of *.exe files


# groups
box_obj  = box.o	boxinit.o	getlen.o	wallloss.o
cond_obj = tmcond.o 	dmdt_int.o  ezcond.o    ezmass.o    \
           getCondSink.o      getMassSink.o     soacond.o
coag_obj = multicoag.o
nuc_obj  = nucleation.o	ricco_nucl.o	dunne_inorg_nucl.o
size_obj = initbounds.o	loginit.o	loginit_bg.o	report.o	\
	   storenm.o	ezwatereqm.o	eznh3eqm.o	\
	   mnfix.o	waterso4.o	waternacl.o  layer_diff.o \
		 waterorg.o
prop_obj = gasdiff.o	aerodens.o
som_obj  = saprc.o	integr2.o	rhs.o	saprc_init.o  saprc14_rev1.o \
					 r_step_autoox.o f_calc_db.o

# executables
box.exe : $(box_obj) $(cond_obj) $(coag_obj) $(nuc_obj) $(size_obj) $(prop_obj) $(som_obj) sizecode.COM modlspc.h gaskin.h common.inc Makefile
	$(comp)  $(exe_flags) -o box.exe $(box_obj) $(cond_obj) $(coag_obj) $(nuc_obj) $(size_obj) $(prop_obj) $(som_obj)

# object files
box.o : box.f sizecode.COM modlspc.h gaskin.h common.inc Makefile
	$(comp) $(o_flags) box.f
wallloss.o: wallloss.f sizecode.COM Makefile
	$(comp) $(o_flags) wallloss.f
boxinit.o : boxinit.f sizecode.COM Makefile
	$(comp) $(o_flags) boxinit.f
tmcond.o : tmcond.f sizecode.COM Makefile
	$(comp) $(o_flags) tmcond.f
ezcond.o : ezcond.f sizecode.COM Makefile
	$(comp) $(o_flags) ezcond.f
ezmass.o : ezmass.f sizecode.COM Makefile
	$(comp) $(o_flags) ezmass.f
getCondSink.o : getCondSink.f sizecode.COM Makefile
	$(comp) $(o_flags) getCondSink.f
getMassSink.o : getMassSink.f sizecode.COM Makefile
	$(comp) $(o_flags) getMassSink.f
getlen.o : getlen.f Makefile
	$(comp) $(o_flags) getlen.f
multicoag.o : multicoag.f sizecode.COM Makefile
	$(comp) $(o_flags) multicoag.f
soacond.o : soacond.f sizecode.COM Makefile
	$(comp) $(o_flags) soacond.f
dmdt_int.o : dmdt_int.f Makefile
	$(comp) $(o_flags) dmdt_int.f
initbounds.o : initbounds.f sizecode.COM Makefile
	$(comp) $(o_flags) initbounds.f
loginit.o : loginit.f sizecode.COM Makefile
	$(comp) $(o_flags) loginit.f
loginit_bg.o : loginit_bg.f sizecode.COM Makefile
	$(comp) $(o_flags) loginit_bg.f
report.o : report.f sizecode.COM Makefile
	$(comp) $(o_flags) report.f
storenm.o : storenm.f sizecode.COM Makefile
	$(comp) $(o_flags) storenm.f
ezwatereqm.o : ezwatereqm.f sizecode.COM Makefile
	$(comp) $(o_flags) ezwatereqm.f
eznh3eqm.o : eznh3eqm.f sizecode.COM Makefile
	$(comp) $(o_flags) eznh3eqm.f
gasdiff.o  : gasdiff.f Makefile
	$(comp) $(o_flags) gasdiff.f
aerodens.o : aerodens.f Makefile
	$(comp) $(o_flags) aerodens.f
mnfix.o : mnfix.f sizecode.COM Makefile
	$(comp) $(o_flags) mnfix.f
waterso4.o : waterso4.f sizecode.COM Makefile
	$(comp) $(o_flags) waterso4.f
waternacl.o : waternacl.f Makefile
	$(comp) $(o_flags) waternacl.f
waterorg.o : waterorg.f sizecode.COM Makefile
	$(comp) $(o_flags) waterorg.f
layer_diff.o : layer_diff.f sizecode.COM Makefile
	$(comp) $(o_flags) layer_diff.f 
nucleation.o : nucleation.f sizecode.COM Makefile
	$(comp) $(o_flags) nucleation.f
ricco_nucl.o  : ricco_nucl.f Makefile
	$(comp) $(o_flags) ricco_nucl.f
dunne_inorg_nucl.o : dunne_inorg_nucl.f Makefile
	$(comp) $(o_flags) dunne_inorg_nucl.f
saprc.o : saprc.f modlspc.h Makefile
	$(comp) $(o_flags) saprc.f
integr2.o : integr2.f modlspc.h common.inc Makefile
	$(comp) $(o_flags) integr2.f
rhs.o : rhs.f modlspc.h gaskin.h common.inc Makefile
	$(comp) $(o_flags) rhs.f
saprc_init.o : saprc_init.f sizecode.COM modlspc.h common.inc gaskin.h Makefile
	$(comp) $(o_flags) saprc_init.f
saprc14_rev1.o : saprc14_rev1.f Makefile
	$(comp) $(o_flags) saprc14_rev1.f
r_step_autoox.o : r_step_autoox.f sizecode.COM modlspc.h common.inc gaskin.h Makefile
	$(comp) $(o_flags) r_step_autoox.f
f_calc_db.o : f_calc_db.f sizecode.COM modlspc.h gaskin.h Makefile 
	$(comp) $(o_flags) f_calc_db.f

# Other Utilities 
# ===============================================================================
clean :
	rm *.o box.exe

clear : 
	rm ../outputs/* 
