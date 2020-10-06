# ----------------------------machine dependent file

# define MEX, in mac mex may also refer to pdftex, a command conflict
MEX = /usr/local/MATLAB/R2019a/bin/mex

# Tested the flag on R2018a
FLAGS = 

FORTRAN   = gfortran
OPTS      = -O3
FPPFLAGS  = -DLAP -DSLAP -fPIC
FOPTS     = -w -Ofast -cpp

LIBNAME = -lfftw3 -lblas 
