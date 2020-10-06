
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
mwSize       :: nt,nc,nw,nints

real*8, allocatable :: t(:),cheb(:),w(:),S(:,:),ab(:,:),st(:,:),x(:)
integer*4, allocatable :: nu(:),nu1(:,:)
integer*4 i,j,k,temp,start,loc,flag,toM,toN,Temp1
real*8 sum1

nt = mxGetM(prhs(1))
nc = mxGetM(prhs(2))
nw = mxGetM(prhs(3))
nints = mxGetN(prhs(4))

allocate(cheb(nc),t(nt),w(nw),ab(2,nints))

call mxCopyPtrToReal8(mxGetPr(prhs(1)),t,nt)
call mxCopyPtrToReal8(mxGetPr(prhs(2)),cheb,nc)
call mxCopyPtrToReal8(mxGetPr(prhs(3)),w,nw)
call mxCopyPtrToReal8(mxGetPr(prhs(4)),ab,nints*2)

allocate(st(nt,nints),nu(nints),nu1(2,nints))
st = 0.0d0
nu = 0
loc = 1
do i = 1,nt
   flag = 0
   start = loc
   do while (flag .eq. 0)
      if (ab(1,start) .lt. t(i) .AND. t(i) .le. ab(2,start)) then
         loc = start
         flag = 1
      else
         start = start + 1
      endif
   end do
   nu(loc) = nu(loc) + 1
   st(nu(loc),loc) = t(i)
end do

nu1 = 0
Temp1 = 0
do i = 1,nints
   if (nu(i) .ne. 0) then
      Temp1 = Temp1 + 1
      nu1(1,Temp1) = nu(i)
      nu1(2,Temp1) = i
   endif
end do

allocate(S(nt,nints*nc),x(nc))
S = 0.0d0
toM = 0
do k = 1,Temp1
x = cheb*(ab(2,nu1(2,k))-ab(1,nu1(2,k)))/2+(ab(1,nu1(2,k))+ab(2,nu1(2,k)))/2
do i = 1,nu1(1,k)
   flag = 0
   do j = 1,nc
      if (abs(st(i,nu1(2,k))-x(j)) .lt. 1.0d-16) then
         temp = j
         flag = 1
         exit
      endif
   enddo
   if (flag .eq. 0) then
      do j = 1,nc
         S(toM+i,(nu1(2,k)-1)*nc+j) = w(j)/(st(i,nu1(2,k))-x(j))
      enddo
      S(toM+i,(nu1(2,k)-1)*nc+1:(nu1(2,k)-1)*nc+nc) = S(toM+i,(nu1(2,k)-1)*nc+1:(nu1(2,k)-1)*nc+nc)  & 
              /sum(S(toM+i,(nu1(2,k)-1)*nc+1:(nu1(2,k)-1)*nc+nc))
   else
      S(toM+i,(nu1(2,k)-1)*nc+1:(nu1(2,k)-1)*nc+nc) = 0.0d0
      S(toM+i,(nu1(2,k)-1)*nc+temp) = 1.0d0
   endif
enddo
toM = toM + nu1(1,k)
enddo
plhs(1) = mxCreateDoubleMatrix(nt, nints*nc, 0)
call mxCopyReal8ToPtr(S, mxGetPr(plhs(1)),nt*nints*nc)
plhs(2) = mxCreateDoubleMatrix(nt, nints, 0)
call mxCopyReal8ToPtr(st, mxGetPr(plhs(2)),nt*nints)
plhs(3) = mxCreateDoubleMatrix(1, 1, 0)
sum1 = sum(nu)
call mxCopyReal8ToPtr(sum1, mxGetPr(plhs(3)),1)
deallocate(x,t,w,S,cheb,ab,st,nu)
end subroutine
