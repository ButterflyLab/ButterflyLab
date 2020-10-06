#include "fintrf.h"

! standard header for mex functions
subroutine mexfunction(nlhs, plhs, nrhs, prhs)
use utils
use chebyshev
implicit double precision (a-k,o-z)

! number of input arguments, number of output arguments
integer      :: nlhs, nrhs
! pointer to inputs and outputs
mwPointer    :: plhs(*), prhs(*)
! get some of the matlab mex functions
mwPointer    :: mxGetPr, mxGetPi, mxCreateDoubleMatrix
! define a size integer so that we can get its type
mwSize       :: n,n1,n2,m1,m2,nts,nnu

type(chebexps_data)           :: chebdata
real*8, allocatable :: twhts(:),ab(:,:),t(:),k(:)
complex*16, allocatable :: r(:),m(:,:)
real*8, allocatable :: psivals(:),avals(:),nu(:),ier(:)
real*8, allocatable :: ts(:),avals0(:),psivals0(:),xs(:),wghts(:)
integer*4, allocatable :: t1(:),k1(:)
integer*4 it,kk,ii,i,jj
real*8 da,db,flag,pi,dd


pi=acos(-1.0d0)

if (n .lt. (2**12)) then
it = 10
else
it = 28
end if

n = mxGetM(prhs(1))
nts = mxGetM(prhs(2))
nnu = mxGetM(prhs(3))

allocate(ts(nts),nu(nnu))
call mxCopyPtrToReal8(mxGetPr(prhs(2)),ts,nts)
call mxCopyPtrToReal8(mxGetPr(prhs(3)),nu,nnu)
call mxCopyPtrToReal8(mxGetPr(prhs(4)),da,1)
call mxCopyPtrToReal8(mxGetPr(prhs(5)),db,1)
call mxCopyPtrToReal8(mxGetPr(prhs(6)),flag,1)

allocate(xs(nts))
!call jacobi_quad_mod(n,da,db,ts,twhts)
xs = int(ts/2/pi*(n+0.0d0))*2*pi/n

kk  = 16
call chebexps(kk,chebdata)

dd      = n*1.0d0
dd      = min(0.10d0,1/dd)
dd      = log(dd)/log(2.0d0)
nints   = ceiling(-dd)+1
nints   = 2*nints
allocate(ab(2,nints))

call jacobi_phase_disc(nints,ab)

allocate(psivals(kk*nints),avals(kk*nints))
!allocate(r(kk*nints),m(kk*nints,nnu))
allocate(r(kk*nints),m(kk*nints,nnu))
!allocate(avals0(nts),psivals0(nts))

do i=1,nnu
    dnu = nu(i)
    call jacobi_phase(chebdata,dnu,da,db,nints,ab,avals,psivals)
    !call jacobi_phase_eval(chebdata,dnu,da,db,nints,ab,avals,psivals,nts,ts,avals0,psivals0)
    if (flag .lt. 0) then
       r = avals*exp(dcmplx(0,1)*(psivals-dnu*ts))
    else
       r = avals*exp(dcmplx(0,1)*(psivals-dnu*xs))
    end if

    m(:,i) = r
end do

allocate(ier(1))
if (nts .eq. int(kk*nints)) then
ier(1) = 1
else
ier(1) = 0
endif
plhs(1)=mxCreateDoubleMatrix(nts, nnu, 1)
plhs(2)=mxCreateDoubleMatrix(1, 1, 0)
call mxCopyComplex16ToPtr(m, mxGetPr(plhs(1)),mxGetPi(plhs(1)),nts*nnu)
call mxCopyReal8ToPtr(ier, mxGetPr(plhs(2)),1)
deallocate(ab,m,r,ts,nu,xs,avals,psivals)

end subroutine

