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
mwSize       :: n,nn1,nn2

type(chebexps_data)           :: chebdata
real*8, allocatable :: c(:),twhts(:),ab(:,:)
complex*16, allocatable :: r(:),rr(:),ier(:),rrr(:)
real*8, allocatable :: psivals(:),avals(:),rd(:)
real*8, allocatable :: ts(:),avals0(:),psivals0(:)
real*8, allocatable :: psival(:,:),aval(:,:),avals1(:),psivals1(:)
integer*4, allocatable :: rd1(:),rd2(:)
integer*4 k
integer*4 it,i,j
real*8 da,db
complex*16 a

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
nn1 = mxGetM(prhs(5))
nn2 = mxGetM(prhs(6))

if (n .lt. 2**12) then
it = 10
else 
it = 28
end if

allocate(c((n-it)**2),r(nn1*nn2),rr(nn2),rrr(nn1),ts(n),twhts(n))
allocate(avals0(n),psivals0(n),avals1(n),psivals1(n),rd1(nn1),rd2(nn2))

call mxCopyPtrToReal8(mxGetPr(prhs(2)),c,(n-it)**2)
call mxCopyPtrToReal8(mxGetPr(prhs(3)),da,1)
call mxCopyPtrToReal8(mxGetPr(prhs(4)),db,1)
allocate(rd(nn1))
call mxCopyPtrToReal8(mxGetPr(prhs(5)),rd,nn1)
rd1 = int(rd+0.01)
deallocate(rd)
allocate(rd(nn2))
call mxCopyPtrToReal8(mxGetPr(prhs(6)),rd,nn2)
rd2 = int(rd+0.01)
deallocate(rd)

call date_and_time(date,time,zone,values1)
call jacobi_quad_mod(n,da,db,ts,twhts)
!ier(1)=1
k  = 16
call chebexps(k,chebdata)

dd      = n*1.0d0
dd      = min(0.10d0,1/dd)
dd      = log(dd)/log(2.0d0)
nints   = ceiling(-dd)+1
nints   = 2*nints
allocate(ab(2,nints))

call jacobi_phase_disc(nints,ab)
!ier(2)=1
allocate(psivals(k*nints),avals(k*nints))
allocate(psival(k*nints,n-it),aval(k*nints,n-it))

do i=it,n-1
dnu = i
call jacobi_phase(chebdata,dnu,da,db,nints,ab,avals,psivals)
psival(:,i-it+1) = psivals
aval(:,i-it+1) = avals
end do
call date_and_time(date,time,zone,values2)
time1=sum((values2(5:8)-values1(5:8))*arr)

r=0
do i=it,n-1
dnu = i
call jacobi_phase_eval(chebdata,dnu,da,db,nints,ab,aval(:,i-it+1),psival(:,i-it+1),n,ts(rd1),avals0,psivals0)
rrr = avals0*exp(dcmplx(0,1)*psivals0)
   do j=it,n-1
      dnu1 = j
      call jacobi_phase_eval(chebdata,dnu1,da,db,nints,ab,aval(:,j-it+1),psival(:,j-it+1),n,ts(rd2),avals1,psivals1)

      rr=avals1*exp(dcmplx(0,1)*psivals1)
      do k=1,nn1
         !r((k-1)*n+1:k*n)=c(rd2(j)-it+1+(rd1(i)-it)*(n-it))*rrr(k)*rr+r((k-1)*n+1:k*n)
         r((k-1)*nn2+1:k*nn2)=c(j-it+1+(i-it)*(n-it))*rrr(k)*rr+r((k-1)*nn2+1:k*nn2)
      end do
      
   end do
end do
!ier(4)=1
!ier=c(1:5)
plhs(1) = mxCreateDoubleMatrix(nn1*nn2, 1, 1)
!plhs(2) = mxCreateDoubleMatrix(5,1,1)
plhs(2) = mxCreateDoubleMatrix(n,1,0)
plhs(3) = mxCreateDoubleMatrix(1,1,0)
call mxCopyComplex16ToPtr(r, mxGetPr(plhs(1)),mxGetPi(plhs(1)),nn1*nn2)
!ier(5)=1
!call mxCopyComplex16ToPtr(ier,mxGetPr(plhs(2)),mxGetPi(plhs(2)),5)
call mxCopyReal8ToPtr(ts,mxGetPr(plhs(2)),n)
call mxCopyReal8ToPtr(time1,mxGetPr(plhs(3)),1)
deallocate(c,twhts,ts,ab,r,psivals,avals,psivals0,avals0,psival,aval,avals1,psivals1)


end subroutine
