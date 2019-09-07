include make.inc

OS = $(shell uname)

.SUFFIXES: .o .a .x .f90 .F90

ALLOBJ = bessel_mod.o alegendre_mod.o utils.o legendre.o

all : clean ${ALLOBJ} fastALegendre_mex.mex fastBessel_mex.mex alegendreTp_mex.mex alegendrePRoot_mex.mex alegendreQRoot_mex.mex alegendreInverse_mex.mex alegendrePhaseDer_mex.mex alegendreJacobi_mex.mex alegendreNRoots_mex.mex legenQuad_mex.mex 

fastALegendre_mex.mex: fastALegendre_mex.F90
	${MEX} ${FLAGS} fastALegendre_mex.F90 $(ALLOBJ)

alegendreNRoots_mex.mex: alegendreNRoots_mex.F90
	${MEX} ${FLAGS} alegendreNRoots_mex.F90 $(ALLOBJ)

alegendreTp_mex.mex: alegendreTp_mex.F90
	${MEX} ${FLAGS} alegendreTp_mex.F90 $(ALLOBJ)

alegendrePRoot_mex.mex: alegendrePRoot_mex.F90
	${MEX} ${FLAGS} alegendrePRoot_mex.F90 $(ALLOBJ)

alegendrePhaseDer_mex.mex: alegendrePhaseDer_mex.F90
	${MEX} ${FLAGS} alegendrePhaseDer_mex.F90 $(ALLOBJ)

alegendreInverse_mex.mex: alegendreInverse_mex.F90
	${MEX} ${FLAGS} alegendreInverse_mex.F90 $(ALLOBJ)

alegendreJacobi_mex.mex: alegendreJacobi_mex.F90
	${MEX} ${FLAGS} alegendreJacobi_mex.F90 $(ALLOBJ)

alegendreQRoot_mex.mex: alegendreQRoot_mex.F90
	${MEX} ${FLAGS} alegendreQRoot_mex.F90 $(ALLOBJ)

fastBessel_mex.mex: fastBessel_mex.F90
	${MEX} ${FLAGS} fastBessel_mex.F90 $(ALLOBJ)

legenQuad_mex.mex: legenQuad_mex.F90
	${MEX} ${FLAGS} legenQuad_mex.F90 $(ALLOBJ)

LINK_MACRO = $< bessel_mod.o alegendre_mod.o  utils.o legendre.o -o $@

clean : 
	rm -f *.a
	rm -f *.o
	rm -f *.x
	rm -f *.mod
	#rm -f *.mexa64
	#rm -f *.mexmaci64

.F90.o : 
	$(FORTRAN) $(OPTS) $(FPPFLAGS) -c $<  -o $@

.f90.o : 
	$(FORTRAN) $(OPTS) $(FPPFLAGS) -c $<  -o $@

.F90.x : 
	$(FORTRAN) $(OPTS) $(FPPFLAGS) -pg $(LINK_MACRO)

.f90.x : 
	$(FORTRAN) $(OPTS) $(LINK_MACRO)
