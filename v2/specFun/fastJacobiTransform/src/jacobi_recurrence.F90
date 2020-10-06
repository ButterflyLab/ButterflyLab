#include "fintrf.h"

! standard header for mex functions
subroutine mexfunction(nlhs, plhs, nrhs, prhs)
use jacobi_transform
implicit double precision (a-j,o-z)

! number of input arguments, number of output arguments
integer      :: nlhs, nrhs  
! pointer to inputs and outputs
mwPointer    :: plhs(*), prhs(*) 
! get some of the matlab mex functions
mwPointer    :: mxGetPr, mxGetPi, mxCreateDoubleMatrix 
! define a size integer so that we can get its type
mwSize       :: n,nn,nts,nnu

real*8, allocatable :: tts(:),vvals0(:,:)
integer(kind = 4)  ::  it,nn1
real*8 da,db,it1

 real*8 arr(4),time1
 character*8 date
 character*10 time
 character*5 zone
 integer*4 values1(8),values2(8)
 
 arr(1)=3600
 arr(2)=60
 arr(3)=1
 arr(4)=0.001
 
n = mxGetM(prhs(1))

allocate(tts(n))

call mxCopyPtrToReal8(mxGetPr(prhs(1)),tts,n)
call mxCopyPtrToReal8(mxGetPr(prhs(2)),da,1)
call mxCopyPtrToReal8(mxGetPr(prhs(3)),db,1)
call mxCopyPtrToReal8(mxGetPr(prhs(4)),it1,1)
it = it1

nn1 = n
it = it
allocate(vvals0(n,it))
call jacobi_recurrence2(nn1,tts,it,da,db,vvals0)


plhs(1) = mxCreateDoubleMatrix(n, it, 0)
call mxCopyReal8ToPtr(vvals0, mxGetPr(plhs(1)),n*it)

deallocate(tts,vvals0)


end subroutine
