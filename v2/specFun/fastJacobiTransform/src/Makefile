include make.inc

OS = $(shell uname)

.SUFFIXES: .o .a .x  .f90 .F90 .f

ALLOBJ = dfft.o dfftpack.o next235.o 
ALLOBJ1 = utils.o amos.o chebyshev.o jacobi_asym.o jacobi_taylor.o jacobi_phase.o jacobi_quad.o
ALLOBJ2 = gspiv.o orthom.o idecomp.o jacobi_exp.o jacobi_transform.o
ALLOBJ3 = nufft1df90.o nufft2df90.o nufft3df90.o

all : ${ALLOBJ} ${ALLOBJ1} ${ALLOBJ2} ${ALLOBJ3} nufft1dIInyumex.mex jacobiexample.mex chebjacex.mex nufft2dIInyumex.mex directjac2.mex extrjac2.mex directjac3.mex \
	extrjac3.mex nufft3dIInyumex.mex directjac1.mex extrjac1.mex getts.mex directinvjac1.mex jacobi_recurrence.mex interpjac1.mex barcycheby.mex jac_exp_extr.mex interpjac2.mex

nufft1dIInyumex.mex: nufft1dIInyumex.F90
	${MEX} ${FLAGS} nufft1dIInyumex.F90 $(ALLOBJ) nufft1df90.o

nufft2dIInyumex.mex: nufft2dIInyumex.F90
	${MEX} ${FLAGS} nufft2dIInyumex.F90 $(ALLOBJ) nufft2df90.o

nufft3dIInyumex.mex: nufft3dIInyumex.F90
	${MEX} ${FLAGS} nufft3dIInyumex.F90 $(ALLOBJ) nufft3df90.o

getts.mex:getts.F90
	${MEX} ${FLAGS} getts.F90 $(ALLOBJ1) $(LIBNAME) 

directjac1.mex:directjac1.F90
	${MEX} ${FLAGS} directjac1.F90 $(ALLOBJ1) $(ALLOBJ2) $(LIBNAME)

directjac2.mex:directjac2.F90
	${MEX} ${FLAGS} directjac2.F90 $(ALLOBJ1) $(ALLOBJ2) $(LIBNAME)

directjac3.mex:directjac3.F90
	${MEX} ${FLAGS} directjac3.F90 $(ALLOBJ1) $(ALLOBJ2) $(LIBNAME)

extrjac1.mex:extrjac1.F90
	${MEX} ${FLAGS} extrjac1.F90 $(ALLOBJ1) $(LIBNAME)

extrjac2.mex:extrjac2.F90
	${MEX} ${FLAGS} extrjac2.F90 $(ALLOBJ1) $(LIBNAME)

extrjac3.mex:extrjac3.F90
	${MEX} ${FLAGS} extrjac3.F90 $(ALLOBJ1) $(LIBNAME)
	
jacobiexample.mex: jacobiexample.F90
	${MEX} ${FLAGS} jacobiexample.F90 $(ALLOBJ1) $(LIBNAME)

chebjacex.mex: chebjacex.F90
	${MEX} ${FLAGS} chebjacex.F90 $(ALLOBJ1) $(ALLOBJ2) $(LIBNAME) 

directinvjac1.mex:directinvjac1.F90
	${MEX} ${FLAGS} directinvjac1.F90 $(ALLOBJ1) $(LIBNAME)

jacobi_recurrence.mex:jacobi_recurrence.F90
	${MEX} ${FLAGS} jacobi_recurrence.F90 $(ALLOBJ1) $(ALLOBJ2) $(LIBNAME)

interpjac1.mex:interpjac1.F90
	${MEX} ${FLAGS} interpjac1.F90 $(ALLOBJ1) $(LIBNAME)

interpjac2.mex:interpjac2.F90
	${MEX} ${FLAGS} interpjac2.F90 $(ALLOBJ1) $(LIBNAME)

barcycheby.mex:barcycheby.F90
	${MEX} ${FLAGS} barcycheby.F90 

jac_exp_extr.mex:jac_exp_extr.F90
	${MEX} ${FLAGS} jac_exp_extr.F90 $(ALLOBJ1) $(ALLOBJ2) $(LIBNAME)
#libdfftpack.a:
#	cd dfftpack && $(MAKE) clean && $(MAKE) && cp libdfftpack.a ..

LINK_MACRO = $< nufft1dIInyumex.o jacobiexample.o chebjacex.o nufft2dIInyumex.o directjac2.o extrjac2.o directjac3.o extrjac3.o nufft3dIInyumex.o directjac1.o extrjac1.o getts.o directinvjac1.o jacobi_recurrence.o interpjac1.o barcycheby.o jac_exp_extr.o interpjac2.o -o $@

clean : 
	rm -f *.a
	rm -f *.o
	rm -f *.x
	rm -f *.mod
	rm -f *.mexa64

.f.o : 
	$(FORTRAN) $(OPTS) $(FPPFLAGS) -c $<  -o $@ $(LIBNAME)

.F90.o : 
	$(FORTRAN) $(OPTS) $(FPPFLAGS) -c $<  -o $@ $(LIBNAME)

.f90.o : 
	$(FORTRAN) $(OPTS) $(FPPFLAGS) -c $<  -o $@ $(LIBNAME)

.f.x : 
	$(FORTRAN) $(OPTS) $(FPPFLAGS) -pg $(LINK_MACRO)

.F90.x : 
	$(FORTRAN) $(OPTS) $(FPPFLAGS) -pg $(LINK_MACRO)

.f90.x : 
	$(FORTRAN) $(OPTS) $(LINK_MACRO)
