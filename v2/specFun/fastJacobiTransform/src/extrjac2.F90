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
mwSize       :: n,n1,n2,m1,m2

type(chebexps_data)           :: chebdata
real*8, allocatable :: twhts(:),ab(:,:),t(:),k(:)
complex*16, allocatable :: r(:),m(:,:),ier(:),rr(:),tt(:)
real*8, allocatable :: psivals(:),avals(:)
real*8, allocatable :: ts(:),avals0(:),psivals0(:),xs(:)
real*8, allocatable :: avals1(:),psivals1(:)
integer*4, allocatable :: t1(:),k1(:)
integer*4 it,kk,ii,i,jj
real*8 da,db,flag,pi


pi=acos(-1.0d0)

!da=-0.50d0
!db=-0.50d0
allocate(ier(5))
n = mxGetM(prhs(1))
if (n .lt. 2**12) then
it = 10
else 
it = 28
end if

m1 = mxGetN(prhs(2))
m2 = mxGetN(prhs(3))
n1 = mxGetM(prhs(2))
n2 = mxGetM(prhs(3))
n1 = n1*m1
n2 = n2*m2
allocate(t(n1),k(n2),t1(n1),k1(n2),m(n1,n2),r(n**2),rr(n),tt(n))
call mxCopyPtrToReal8(mxGetPr(prhs(2)),t,n1)
call mxCopyPtrToReal8(mxGetPr(prhs(3)),k,n2)
call mxCopyPtrToReal8(mxGetPr(prhs(4)),flag,1)
call mxCopyPtrToReal8(mxGetPr(prhs(5)),da,1)
call mxCopyPtrToReal8(mxGetPr(prhs(6)),db,1)
t1=int(t+0.5)
k1=int(k+0.5)
ier(1)=flag
allocate(ts(n),twhts(n),xs(n))
allocate(avals0(n),psivals0(n),avals1(n),psivals1(n))
ier(1:5)=k1(1:5)
!ier(3:5)=k1(1:3)
call jacobi_quad_mod(n,da,db,ts,twhts)
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
!ier(1)=nints
!ier(2)=ab(1,1)
!ier(3)=ab(2,nints)
!ier(4:5)=0
allocate(psivals(kk*nints),avals(kk*nints))
!ier(2)=1
!ier=xs(1:5)
m=0
do i=1,n2
if (k1(i) .gt. n*it) then
dnu = mod(k1(i),n)
if (abs(dnu) .le. 0.10d0) then
   dnu = n
end if
if (dnu .gt. it) then
    dnu = dnu-1
    
    j = int((k1(i)-0.5)/n)
    call jacobi_phase(chebdata,j,da,db,nints,ab,avals,psivals)
    !ier(1)=j
    !ier(2:5)=psivals(1:4)
    call jacobi_phase_eval(chebdata,j,da,db,nints,ab,avals,psivals,n,ts,avals0,psivals0)
    !ier=avals0(1:5)
    call jacobi_phase(chebdata,dnu,da,db,nints,ab,avals,psivals)
    call jacobi_phase_eval(chebdata,dnu,da,db,nints,ab,avals,psivals,n,ts,avals1,psivals1)
    if (flag .lt. 0) then
    tt = avals0*exp(dcmplx(0,1)*(psivals0-j*ts))
    rr = avals1*exp(dcmplx(0,1)*(psivals1-dnu*ts))
    else
    tt = avals0*exp(dcmplx(0,1)*(psivals0-j*xs))
    rr = avals1*exp(dcmplx(0,1)*(psivals1-dnu*xs))
    end if
    do ii=1,n1
       jj = int(t1(ii)/n)+1
       kk = mod(t1(ii),n)
       if (abs(real(kk)) .lt. 0.1d0) then
          kk = n
          jj = t1(ii)/n
       end if
       m(ii,i) = tt(jj)*rr(kk)
       !r((ii-1)*n+1:ii*n) = tt(ii)*rr
       !r((a-1)*n+1:a*n)=avals0(a)*exp(dcmplx(0,1)*(psivals0(a)-j*ts(a)-j*flag*(xs(a)-ts(a))))*avals1*exp(dcmplx(0,1)   &
        !*(psivals1-dnu*ts-dnu*flag*(xs-ts)))+r((a-1)*n+1:a*n)
    end do
    
    !ier=avals0(1:5)*exp(dcmplx(0,1)*(psivals0(1:5)-j*ts(1:5)-j*flag*(xs(1:5)-ts(1:5))))
!    ier=avals1(1:5)*exp(dcmplx(0,1)*(psivals1(1:5)-dnu*ts(1:5)-dnu*flag*(xs(1:5)-ts(1:5))))
    !m(:,i)=r(t1)
    !ier=r([2,5,6,7,9])
end if
end if
end do

plhs(1)=mxCreateDoubleMatrix(n1, n2, 1)
plhs(2)=mxCreateDoubleMatrix(5,1,1)
call mxCopyComplex16ToPtr(m, mxGetPr(plhs(1)),mxGetPi(plhs(1)),n1*n2)
call mxCopyComplex16ToPtr(ier,mxGetPr(plhs(2)),mxGetPi(plhs(2)),5)
deallocate(twhts,ab,t,k,t1,k1,m,r,ts,xs,avals,psivals,avals0,psivals0,avals1,psivals1)

end subroutine
