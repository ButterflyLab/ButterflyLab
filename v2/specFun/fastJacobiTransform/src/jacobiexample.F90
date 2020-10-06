#include "fintrf.h"

! standard header for mex functions
subroutine mexfunction(nlhs, plhs, nrhs, prhs)
use utils
use chebyshev

implicit double precision (a-h,o-z)

! number of input arguments, number of output arguments
integer      :: nlhs, nrhs  
! pointer to inputs and outputs
mwPointer    :: plhs(*), prhs(*) 
! get some of the matlab mex functions
mwPointer    :: mxGetPr, mxGetPi, mxCreateDoubleMatrix 
! define a size integer so that we can get its type
mwSize       :: nj,nts

type(chebexps_data)           :: chebdata
real*8, allocatable :: ab(:,:)
real*8, allocatable :: psivals(:),avals(:)
real*8, allocatable :: ts(:),avals0(:),psivals0(:),polvals(:),polvals0(:)
real*8, allocatable :: xs(:),twhts(:),tss(:,:)
complex*16,allocatable :: jacobi1(:,:),jacobi2(:,:)
real*8 :: nts1,da,db
integer*4  k,nt
real*8 pi
integer nints,it

!allocate(nts1(2,1))
pi  = acos(-1.0d0)
!call mxCopyPtrToReal8(mxGetPr(prhs(1)),nts1,1)
!nts = ceiling(nts1(1,1))
nts = mxGetM(prhs(1))
if (nts .lt. 2**(12)) then
it = 10
else
it = 28
end if
nt = nts
!nts = 2**9
allocate(jacobi1(nts,nts-it),jacobi2(nts,nts-it))
allocate(ts(nts),xs(nts),twhts(nts),avals0(nts),tss(nts,1))
allocate(psivals0(nts),polvals(nts),polvals0(nts))

!da=0.25d0
!db=0.25d0
call mxCopyPtrToReal8(mxGetPr(prhs(2)),da,1)
call mxCopyPtrToReal8(mxGetPr(prhs(3)),db,1)
call jacobi_quad_mod(nts,da,db,ts,twhts)

xs = int(ts/2/pi*(nts+0.0d0))*2*pi/nts
!call quicksort(nt,ts)
!call quicksort(nt,xs)




k  = 16
call chebexps(k,chebdata)

dd      = nts*1.0d0
dd      = min(0.10d0,1/dd)
dd      = log(dd)/log(2.0d0)
nints   = ceiling(-dd)+1
nints   = 2*nints
allocate(ab(2,nints))

call jacobi_phase_disc(nints,ab)

allocate(psivals(k*nints),avals(k*nints))

do i=it,nts-1


! this must be integer because we use the 
!recurrence relations to test accuracy

dnu = i
!da  = 0.25d0
!db  = 0.25d0
!p   = dnu + (da+db+1)/2
!allocate(xx(0:nu))


!call prini("nints = ",nints)
!call prin2("before jacobi_phase, ab = ",ab)



call jacobi_phase(chebdata,dnu,da,db,nints,ab,avals,psivals)
!call prin2("time to construct phase = ",(t2-t1))

call jacobi_phase_eval(chebdata,dnu,da,db,nints,ab,avals,psivals,nts,ts,avals0,psivals0)
!call prin2("average eval time = ",(t2-t1)/nts)
!polvals0 = cos(psivals0)*avals0
!print *,size(avals0),size(psivals0)
jacobi1(:,i-it+1)=avals0*exp(dcmplx(0,1)*(psivals0-i*ts))
jacobi2(:,i-it+1)=avals0*exp(dcmplx(0,1)*(psivals0-i*xs))
end do

plhs(1) = mxCreateDoubleMatrix(nts, 1, 0)
plhs(2) = mxCreateDoubleMatrix(nts, nts-it, 1)
plhs(3) = mxCreateDoubleMatrix(nts, nts-it, 1)
!tss(:,1) = ts
call mxCopyReal8ToPtr(ts, mxGetPr(plhs(1)),nts) 
call mxCopyComplex16ToPtr(jacobi1, mxGetPr(plhs(2)),mxGetPi(plhs(2)),nts*(nts-it)) 
call mxCopyComplex16ToPtr(jacobi2, mxGetPr(plhs(3)),mxGetPi(plhs(3)),nts*(nts-it)) 

deallocate(jacobi1,jacobi2,ts,xs,twhts,avals0)
deallocate(psivals0,polvals,polvals0)
end subroutine mexfunction
