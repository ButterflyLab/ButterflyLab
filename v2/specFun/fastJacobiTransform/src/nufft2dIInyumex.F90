#include "fintrf.h"

! standard header for mex functions
subroutine mexfunction(nlhs, plhs, nrhs, prhs)

  implicit double precision (a-h,o-z)

  ! number of input arguments, number of output arguments
  integer      :: nlhs, nrhs  
  ! pointer to inputs and outputs
  mwPointer    :: plhs(*), prhs(*) 
  ! get some of the matlab mex functions
  mwPointer    :: mxGetPr, mxGetPi, mxCreateDoubleMatrix 
  ! define a size integer so that we can get its type
  mwSize       :: nj

  ! input variable
  real*8 eps
  integer iflag,ms
  real*8, allocatable ::  xj(:),yj(:)
  complex*16, allocatable ::  fk(:,:)

  ! output variable
  complex*16, allocatable ::  cj(:)

  ! aux variable
  integer ier,iflag1,ms1,n

  nj     = mxGetM(prhs(1))
  n = nj
  n = int(sqrt(real(n)))

  ! copy the right-hand side argument in matlab to dnu and t
  allocate(cj(nj),xj(nj),yj(nj))

  call mxCopyPtrToReal8(mxGetPr(prhs(1)), xj, nj)
  call mxCopyPtrToReal8(mxGetPr(prhs(2)), yj, nj)
  call mxCopyPtrToInteger4(mxGetPr(prhs(3)),iflag,1)
  call mxCopyPtrToReal8(mxGetPr(prhs(4)), eps, 1)
  !call mxCopyPtrToInteger4(mxGetPr(prhs(5)),ms,1)
  !ms1=int(ms)
  allocate(fk(n,n))
  call mxCopyPtrToComplex16(mxGetPr(prhs(5)),mxGetPi(prhs(5)), fk, n*n)
  !iflag1=int(iflag)
  
  !print *,'iflag=',iflag
  !print *,'ms   =',ms
  ! initialize the left-hand side of the function
  plhs(1) = mxCreateDoubleMatrix(nj, 1, 1)


  call nufft2d2f90(nj,xj,yj,cj,1,eps,n,n,fk,ier)


  ! copy the output to the left-hand side
  call mxCopyComplex16ToPtr(cj, mxGetPr(plhs(1)),mxGetPi(plhs(1)),nj) 
  !call mxCopyReal8toPtr(iflag,mxGetPr(plhs(2)),1)
  deallocate(cj,xj,fk,yj)

end subroutine mexfunction
