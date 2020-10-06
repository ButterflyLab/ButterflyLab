subroutine jacobi_solve(k,xscheb,S1,S2,S3,aintl,nints,ab,dnu,da,db,psivals,avals)
use chebyshev
implicit double precision (a-h,o-z)

integer                        :: nints
double precision               :: ab(2,nints),avals(k,nints),psivals(k,nints)
double precision               :: xscheb(k),S1(k,k),S2(k,k),S3(k,k),aintl(k,k)

!
!  Construct the nonoscillatory phase and amplitude functions for Jacobi's differential
!  equation.  These functions are represented via piecewise Chebyshev expansions given
!  on a collection of intervals.
!
!  Input parameters:
!
!  Output parameters:
!

double precision      :: valsmp(k,nints),valsmpp(k,nints),alphap(k,nints)
data  pi        / 3.14159265358979323846264338327950288d0 /
data  piovertwo / 1.57079632679489661923132169163975144d0 /

p = dnu + (da+db+1)/2

!
!  Compute the value of the amplitude function and its first two derivatives
!  at the point pi/2.
!

call jacobi_asym3(dnu,da,db,c3,c2,c1)


!
!  Solve the third order linear ODE for the square of the amplitude function
!

do int=nints,1,-1
a = ab(1,int)
b = ab(2,int)
call jacobi_solve_int(dnu,da,db,a,b,k,xscheb,S1,S2,S3, &
  c1,c2,c3,avals(:,int),valsmp(:,int),valsmpp(:,int))
c1 = valsmpp(1,int)
c2 = valsmp(1,int)
c3 = avals(1,int)
end do


!
!  Record the values of alpha' and alpha'' at the left endpoint for later
!  user
!

apval  = 1/avals(1,1)
appval = -valsmp(1,1)/avals(1,1)**2

!
!  Construct the phase function and translate it appropriately
!

alphap = 1.0d0/avals
avals  = sqrt(avals)


c1     = 0
do int=1,nints
a              = ab(1,int)
b              = ab(2,int)
psivals(:,int) = c1 + matmul((b-a)/2*aintl,alphap(:,int))
c1             = psivals(k,int)
end do


a = ab(1,1)
call jacobi_asym0(dnu,da,db,a,y,yp)
c1 = y*sqrt(apval)
c2 = yp/sqrt(apval)+y*appval/(2*apval**(1.5d0))
a2  = atan2(c1,c2)

if(a2 .gt. pi) then
a2 = a2-pi
endif

if(a2 .le. 0) then
a2 = a2 + pi
endif

psivals = psivals+a2-piovertwo

!
!  Subtract p*t
!

do int=1,nints
a = ab(1,int)
b = ab(2,int)
do i=1,k
t = (b-a)/2 * xscheb(i) + (b+a)/2
psivals(i,int) = psivals(i,int)   - p*t
end do
end do
 

psivals = psivals/p


end subroutine



subroutine jacobi_solve_int(dnu,dalpha,dbeta,a,b,k,xs,S1,S2,S3,c1,c2,c3,valsm,valsmp,valsmpp)
use utils
use chebyshev
implicit double precision (a-h,o-z)

double precision       :: S1(k,k),S2(k,k),S3(k,k),xs(k)
double precision       :: valsm(k),valsmp(k),valsmpp(k)

double precision :: qvals(k),qpvals(k),amatr(k,k),xx(k),xx2(k)
double precision :: sigma(k)

dd   = (b-a)/2
xx = xs*(b-a)/2 + (a-b)/2
xx2 = 0.125d0 * (a-b)**2 * (xs-1)**2

!
!  Evaluate Q and its derivative at the specified points
!

do i=1,k
t         = (b-a)/2 * xs(i) + (b+a)/2
qvals(i)  = dnu*(1+dalpha+dbeta+dnu)+((1+dalpha+dbeta-1.d0*(-dalpha+  &
dbeta)*Cos(t))*1/sin(t)**2)/2.d0-((1+dalpha+dbeta)*cos(t)/sin(t)-1.d0*(-dalpha+dbeta)*&
  1/sin(t))**2/4.d0

qpvals(i) = (((-2+4.d0*dalpha**2+4.d0*dbeta**2)*Cos(t)+(dalpha-1.d0*dbeta)*(dalpha  &
+dbeta)*(3+Cos(2.d0*t)))*1/sin(t)**3)/4.d0
end do

!
!  Form the linear system to solve
!

do i=1,k
amatr(i,:) = 4*qvals(i) *S2(i,:)*dd**2 + 2 * qpvals(i) * S3(i,:)*dd**3
end do
do i=1,k
amatr(i,i) = amatr(i,i)+ 1.0d0
end do

sigma = -4*qvals * ( c2 + c1 * xx ) - 2 * qpvals * (c3 + c2 * xx + c1 * xx2)

!
!  Solve it 
!

call qrsolv(amatr,k,sigma)

valsmpp = c1 + matmul(S1*dd,sigma)
valsmp  = c2 + matmul(S1*dd,valsmpp)
valsm   = c3 + matmul(S1*dd,valsmp)

end subroutine





! program test_solve_jacobi
! use utils
! use chebyshev
! implicit double precision (a-h,o-z)

! type(chebexps_data)           :: chebdata
! double precision, allocatable :: avals(:,:),psivals(:,:)
! double precision, allocatable :: ab(:,:),valsm(:,:),valsmp(:,:),valsmpp(:,:)
! double precision, allocatable :: S1(:,:),S2(:,:),S3(:,:),coefs(:,:)

! dnu    = 100.0d0
! da     = 0.25d0
! db     = 0.00d0
! pi     = acos(-1.0d0)
! k      = 30

! call chebexps2(k,chebdata)

! allocate(S1(k,k),S2(k,k),S3(k,k))
! S1 = chebdata%aintr
! S2 = matmul(S1,S1)
! S3 = matmul(S1,S2)

! nints = 13
! allocate(ab(2,nints),valsm(k,nints),valsmp(k,nints),valsmpp(k,nints),coefs(k,nints))
! allocate(avals(k,nints),psivals(k,nints))


! !
! !  Set the list of intervals
! !
! p = dnu+(da+db+1)/2
! a = 2/p
! b = pi/2

! do int=1,nints-1
! ab(1,nints-1-int+2) = a + (b-a) * 2.0d0**(-int) 
! ab(2,nints-1-int+2) = a + (b-a) * 2.0d0**(1-int) 
! end do

! ab(1,1) = a
! ab(2,1) = ab(1,2)

! call elapsed(t1)
! call jacobi_solve(chebdata%k,chebdata%xs,S1,S2,S3,chebdata%aintl,nints,ab,dnu,da,db,psivals,avals)
! call elapsed(t2)
! print *,"time = ", (t2-t1)
! print *,""

! x = 1.0d0
! call chebpw_eval(nints,ab,k,chebdata%xs,avals,x,aval)
! call chebpw_eval(nints,ab,k,chebdata%xs,psivals,x,psival)

! print *,aval,0.0996882556288365776669641708501876308d0
! print *,psival
! print *,aval*cos(psival)
! print *,0.0465831293979657797330681514249250580d0
! stop

! ! call jacobi_asym3(dnu,da,db,c3,c2,c1)
! ! do int=nints,1,-1
! ! a = ab(1,int)
! ! b = ab(2,int)
! ! call jacobi_solve_int(dnu,da,db,a,b,k,chebdata%xs,S1,S2,S3,&
! !   c1,c2,c3,valsm(:,int),valsmp(:,int),valsmpp(:,int))
! ! call elapsed(t2)
! ! c1 = valsmpp(1,int)
! ! c2 = valsmp(1,int)
! ! c3 = valsm(1,int)
! ! end do
! ! print *,"time = ",t2-t1

! ! print *,c1
! ! print *,c2
! ! print *,c3

! end program
