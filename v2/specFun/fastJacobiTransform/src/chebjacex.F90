#include "fintrf.h"

! standard header for mex functions
subroutine mexfunction(nlhs, plhs, nrhs, prhs)

use utils
use chebyshev
use idecomp
use jacobi_exp
use jacobi_transform

  implicit double precision (a-h,o-z)

  ! number of input arguments, number of output arguments
  integer      :: nlhs, nrhs  
  ! pointer to inputs and outputs
  mwPointer    :: plhs(*), prhs(*) 
  ! get some of the matlab mex functions
  mwPointer    :: mxGetPr, mxGetPi, mxCreateDoubleMatrix 
  ! define a size integer so that we can get its type
  mwSize       :: dmax


type(jacobi_expansion_data)   :: expdata
type(jacobi_transform_data)   :: jacdata
real*8,allocatable :: r(:,:),ts(:)
!integer*4,allocatable :: idxs(:)
complex*16,allocatable :: expvals(:,:)
real*8   da,db,eps,rdmax
integer*4  n,iffactor
integer krank


iffactor=1
dmax = mxGetM(prhs(1))
n = dmax
rdmax = dmax
call mxCopyPtrToReal8(mxGetPr(prhs(2)),da,1)
call mxCopyPtrToReal8(mxGetPr(prhs(3)),db,1)
call mxCopyPtrToReal8(mxGetPr(prhs(4)),eps,1)

call jacobi_expansion(eps,iffactor,rdmax,da,db,expdata)
call jacobi_transform_prepare(expdata,n,jacdata)
krank = jacdata%krank
allocate(r(dmax,krank),expvals(dmax,krank),ts(dmax))
r = jacdata%r
ts = jacdata%ts
expvals = jacdata%expvals

!plhs(1) = mxCreateDoubleMatrix(1,1,0)
plhs(1) = mxCreateDoubleMatrix(dmax, krank, 0)
plhs(2) = mxCreateDoubleMatrix(dmax, krank, 1)
plhs(3) = mxCreateDoubleMatrix(dmax,1,0)

!call mxCopyInteger4ToPtr(krank,mxGetPr(plhs(1)),1)
call mxCopyReal8ToPtr(r,mxGetPr(plhs(1)),dmax*krank)
call mxCopyComplex16ToPtr(expvals,mxGetPr(plhs(2)),mxGetPi(plhs(2)),dmax*krank)
call mxCopyReal8ToPtr(ts,mxGetPr(plhs(3)),dmax)
deallocate(r,expvals,ts)
end subroutine mexfunction
