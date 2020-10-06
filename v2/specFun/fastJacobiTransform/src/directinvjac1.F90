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
mwSize       :: n,nn,nts,nnu

type(chebexps_data)           :: chebdata
real*8, allocatable :: c(:),twhts(:),ab(:,:)
complex*16, allocatable :: ier(:)
real*8, allocatable :: psivals(:),avals(:),nu(:),r(:)
real*8, allocatable :: ts(:),avals0(:),psivals0(:),wghts(:)
real*8, allocatable :: psival(:,:),aval(:,:),rd(:)
integer*4, allocatable :: rd1(:)
integer*4 k,ii,jj,kk
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
nn = mxGetM(prhs(5))
nts = mxGetM(prhs(6))
nnu = mxGetM(prhs(7))

if (n .lt. (2**12)) then
it = 10
else 
it = 28
end if

allocate(c(nts),r(nn),ts(nts),nu(nnu),twhts(n),rd(nn),rd1(nn),wghts(nts))
allocate(avals0(nts),psivals0(nts))

call mxCopyPtrToReal8(mxGetPr(prhs(2)),c,nts)
call mxCopyPtrToReal8(mxGetPr(prhs(3)),da,1)
call mxCopyPtrToReal8(mxGetPr(prhs(4)),db,1)
call mxCopyPtrToReal8(mxGetPr(prhs(5)),rd,nn)
rd1 = int(rd+0.01)
call mxCopyPtrToReal8(mxGetPr(prhs(6)),ts,nts)
call mxCopyPtrToReal8(mxGetPr(prhs(7)),nu,nnu)
call mxCopyPtrToReal8(mxGetPr(prhs(8)),wghts,nts)

call date_and_time(date,time,zone,values1)
!call jacobi_quad_mod(n,da,db,ts,twhts)
!ier(1)=1
k  = 16
call chebexps(k,chebdata)

dd      = n*1.0d0
dd      = min(0.10d0,1/dd)
dd      = log(dd)/log(2.0d0)
nints   = ceiling(-dd)+1
nints   = 2*nints
allocate(ab(2,nints))
wghts = sqrt(wghts)

call jacobi_phase_disc(nints,ab)
!ier(2)=1
allocate(psivals(k*nints),avals(k*nints))
allocate(psival(k*nints,nn),aval(k*nints,nn))

do i=1,nn
dnu = rd(i)+it-1
call jacobi_phase(chebdata,dnu,da,db,nints,ab,avals,psivals)
psival(:,i) = psivals
aval(:,i) = avals
end do

call date_and_time(date,time,zone,values2)
time1=sum((values2(5:8)-values1(5:8))*arr)



!r=0
do i=1,nn
dnu = rd(i)+it-1
call jacobi_phase_eval(chebdata,dnu,da,db,nints,ab,aval(:,i),psival(:,i),nts,ts,avals0,psivals0)
r(i) = sum(avals0*dcos(psivals0)*c*wghts)
end do




plhs(1) = mxCreateDoubleMatrix(nn, 1, 0)
!plhs(2) = mxCreateDoubleMatrix(5,1,1)
!plhs(2) = mxCreateDoubleMatrix(n,1,0)
plhs(2) = mxCreateDoubleMatrix(1,1,0)
call mxCopyReal8ToPtr(r, mxGetPr(plhs(1)),nn)
!ier(5)=1
!call mxCopyComplex16ToPtr(ier,mxGetPr(plhs(2)),mxGetPi(plhs(2)),5)
!call mxCopyReal8ToPtr(ts,mxGetPr(plhs(2)),n)
call mxCopyReal8ToPtr(time1,mxGetPr(plhs(2)),1)
deallocate(c,twhts,ts,ab,r,psivals,avals,psivals0,avals0,psival,aval)


end subroutine
