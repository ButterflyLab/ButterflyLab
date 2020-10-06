
#include "fintrf.h"

! standard header for mex functions
subroutine mexfunction(nlhs, plhs, nrhs, prhs)
use utils
use chebyshev
implicit double precision (a-j,o-z)

! number of input arguments, number of output arguments
integer      :: nlhs, nrhs  
! pointer to inputs and outputs
mwPointer    :: plhs(*), prhs(*) 
! get some of the matlab mex functions
mwPointer    :: mxGetPr, mxGetPi, mxCreateDoubleMatrix 
! define a size integer so that we can get its type
mwSize       :: n


real*8, allocatable :: twhts(:)
real*8, allocatable :: ts(:)
real*8 da,db

n = mxGetM(prhs(1))


allocate(ts(n),twhts(n))


call mxCopyPtrToReal8(mxGetPr(prhs(2)),da,1)
call mxCopyPtrToReal8(mxGetPr(prhs(3)),db,1)


call jacobi_quad_mod(n,da,db,ts,twhts)

plhs(1) = mxCreateDoubleMatrix(n, 1, 0)
plhs(2) = mxCreateDoubleMatrix(n, 1, 0)

call mxCopyReal8ToPtr(ts, mxGetPr(plhs(1)),n)
call mxCopyReal8ToPtr(twhts, mxGetPr(plhs(2)),n)

deallocate(twhts,ts)


end subroutine
