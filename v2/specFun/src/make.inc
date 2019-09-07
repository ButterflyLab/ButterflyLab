# ----------------------------machine dependent file

# define MEX, in mac mex may also refer to pdftex, a command conflict
MEX = /Applications/MATLAB_R2018b.app/bin/mex
#MEX = /usr/local/MATLAB/R2017a/bin/mex
#MEX = mex

# Tested the flag on R2017a
FLAGS = -compatibleArrayDims -g “COMPFLAGS=$COMPFLAGS-fdefault-integer-8”

FORTRAN   = gfortran
OPTS      = -O3
FPPFLAGS  = -DLAP -DSLAP -fPIC

LIBNAME = specFunc
