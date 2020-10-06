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
complex*16, allocatable :: r(:),m(:,:),ier(:),rr(:),tt(:),ee(:)
real*8, allocatable :: psivals(:),avals(:),avals2(:),psivals2(:)
real*8, allocatable :: ts(:),avals0(:),psivals0(:),xs(:)
real*8, allocatable :: avals1(:),psivals1(:)
integer*4, allocatable :: t1(:),k1(:)
integer*4 it,kk,ii,i,jj,mm,nn
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
allocate(t(n1),k(n2),t1(n1),k1(n2),m(n1,n2),r(n**3),rr(n),tt(n),ee(n))
call mxCopyPtrToReal8(mxGetPr(prhs(2)),t,n1)
call mxCopyPtrToReal8(mxGetPr(prhs(3)),k,n2)
call mxCopyPtrToReal8(mxGetPr(prhs(4)),flag,1)
call mxCopyPtrToReal8(mxGetPr(prhs(5)),da,1)
call mxCopyPtrToReal8(mxGetPr(prhs(6)),db,1)

t1=int(t+0.01)
k1=int(k+0.01)
ier(1)=flag
allocate(ts(n),twhts(n),xs(n))
allocate(avals0(n),psivals0(n),avals1(n),psivals1(n),avals2(n),psivals2(n))
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
if (k1(i) .gt. it*(n**2)) then
dnu = mod(k1(i),n**2)
if (abs(dnu) .le. 0.10d0) then
   dnu = n**2
end if
if (dnu .gt. it*n) then
  dnu1 = mod(int(dnu),n)
  if (abs(dnu1) .le. 0.10d0) then
   dnu1 = n
  end if
  if (dnu1 .gt. it) then
    !dnu1 = dnu1-1
   
    a = real(int((k1(i)-0.5)/(n**2)))
    call jacobi_phase(chebdata,a,da,db,nints,ab,avals,psivals)
    call jacobi_phase_eval(chebdata,a,da,db,nints,ab,avals,psivals,n,ts,avals0,psivals0)
    b = real(int((dnu-0.5)/n))
    call jacobi_phase(chebdata,b,da,db,nints,ab,avals,psivals)
    call jacobi_phase_eval(chebdata,b,da,db,nints,ab,avals,psivals,n,ts,avals1,psivals1)
    dnu1 = dnu1-1
    call jacobi_phase(chebdata,dnu1,da,db,nints,ab,avals,psivals)
    call jacobi_phase_eval(chebdata,dnu1,da,db,nints,ab,avals,psivals,n,ts,avals2,psivals2)
    if (flag .lt. 0) then
       tt = avals0*exp(dcmplx(0,1)*(psivals0-a*ts))
       rr = avals1*exp(dcmplx(0,1)*(psivals1-b*ts))
       ee = avals2*exp(dcmplx(0,1)*(psivals2-dnu1*ts))
    else
       tt = avals0*exp(dcmplx(0,1)*(psivals0-a*xs))
       rr = avals1*exp(dcmplx(0,1)*(psivals1-b*xs))
       ee = avals2*exp(dcmplx(0,1)*(psivals2-dnu1*xs))
    end if
    do ii=1,n1
       jj = int(t1(ii)/(n**2))+1
       kk = mod(t1(ii),n**2)
       if (abs(real(kk)) .lt. 0.10d0) then
          jj = t1(ii)/(n**2)
          kk = n
          mm = n
       else
          nn = t1(ii)-(jj-1)*(n**2)
          kk = int(nn/n)+1
          mm = mod(nn,n)
          if (abs(real(mm)) .lt. 0.10d0) then
             kk = nn/n
             mm = n
          else
             mm = nn - (kk-1)*n
          end if
       end if
       m(ii,i) = tt(jj)*rr(kk)*ee(mm)
    end do
    
  end if
end if
end if
end do

plhs(1)=mxCreateDoubleMatrix(n1, n2, 1)
plhs(2)=mxCreateDoubleMatrix(5,1,1)
call mxCopyComplex16ToPtr(m, mxGetPr(plhs(1)),mxGetPi(plhs(1)),n1*n2)
call mxCopyComplex16ToPtr(ier,mxGetPr(plhs(2)),mxGetPi(plhs(2)),5)
deallocate(twhts,ab,t,k,t1,k1,m,r,ts,xs,avals,psivals,avals0,psivals0,avals1,psivals1,tt,rr,ee)
deallocate(avals2,psivals2)

end subroutine
