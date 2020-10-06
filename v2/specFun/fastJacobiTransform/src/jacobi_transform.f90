!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  This module contains code for applying the forward and backward Jacobi transforms
!  rapidly.  Before transforms can be applied
!
!  For our purposes, the forward Jacobi transform is the mapping which takes the
!  coefficients in an expansion to 
!
!
!
!
!  The following routines should be regarded as publicly callable:
!
!   jacobi_transform_prepare - 
!
!   jacobi_transform_size - return the size of the 
!
!   jacobi_transform_destroy - free the memory
!
!   jacobi_transform_forward - 
!
!   jacobi_transform_backward - 
!
!   jacobi_transform_forward_bf - apply the foward transform by brute force
!
!   jacobi_transform_backward_bf - apply the backward transform by brute force
!
!   jacobi_transform_write - write 
!
!   jacobi_transform_read -   
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


#define USE_FFTW

module jacobi_transform

use utils
use chebyshev
use jacobi_exp


type      jacobi_transform_data 
integer                         :: n,krank,nrec

double precision                :: da,db
integer*8                       :: iplan
double precision, allocatable   :: r(:,:),vals0(:,:),ts(:),whts(:),dnus(:),xs(:),dnus0(:)
integer, allocatable            :: idxs(:)
double complex, allocatable     :: expvals(:,:)
double complex, pointer         :: z1(:),z2(:),z3(:)
double complex, allocatable     :: wsave(:)
end type  jacobi_transform_data

double precision, private       :: pi
data pi            / 3.14159265358979323846264338327950288d0 /


INTEGER, private :: FFTW_FORWARD
PARAMETER (FFTW_FORWARD=-1)
INTEGER, private :: FFTW_BACKWARD
PARAMETER (FFTW_BACKWARD=+1)
INTEGER, private :: FFTW_MEASURE
PARAMETER (FFTW_MEASURE=0)
INTEGER, private :: FFTW_ESTIMATE
PARAMETER (FFTW_ESTIMATE=64)


contains


subroutine jacobi_transform_destroy(jacdata)
implicit double precision (a-h,o-z)
type(jacobi_transform_data)     :: jacdata

deallocate(jacdata%r,jacdata%vals0,jacdata%ts,jacdata%whts,jacdata%dnus,jacdata%xs)
deallocate(jacdata%dnus0,jacdata%idxs,jacdata%expvals,jacdata%z1,jacdata%z2)
deallocate(jacdata%z3)

#ifdef USE_FFTW
call dfftw_destroy_plan(jacdata%iplan,n,jacdata%z1,jacdata%z2,FFTW_BACKWARD,FFTW_ESTIMATE)
#else
deallocate(jacdata%wsave)
#endif

end subroutine

subroutine jacobi_transform_size(jacdata,dsize)
implicit double precision (a-h,o-z)
type(jacobi_transform_data)     :: jacdata

dsize = 12 + 8 
dsize = dsize + 8  * ( size(jacdata%r) )
dsize = dsize + 8  * ( size(jacdata%vals0) )
dsize = dsize + 8  * ( size(jacdata%ts) )
dsize = dsize + 8  * ( size(jacdata%whts) )
dsize = dsize + 8  * ( size(jacdata%dnus) )
dsize = dsize + 8  * ( size(jacdata%xs) )
dsize = dsize + 8  * ( size(jacdata%dnus0) )
dsize = dsize + 4  * ( size(jacdata%idxs) )
dsize = dsize + 16 * ( size(jacdata%expvals) )
dsize = dsize + 16 * ( size(jacdata%z1) )
dsize = dsize + 16 * ( size(jacdata%z2) )
dsize = dsize + 16 * ( size(jacdata%z3) )

#ifndef USE_FFTW
dsize = dsize + 16 * ( size(jacdata%wsave) )
#endif

dsize = dsize/(1024.0d0*1024.0d0)

! print *,""
! print *,size(jacdata%r)*8/(1024d0*1024d0)
! print *,size(jacdata%vals0)*8/(1024d0*1024d0)
! print *,size(jacdata%expvals)*16/(1024d0*1024d0)
! print *,size(jacdata%z1)*16/(1024d0*1024d0)
! print *,size(jacdata%z2)*16/(1024d0*1024d0)
! print *,size(jacdata%z3)*16/(1024d0*1024d0)
! print *,""

! print *,dsize
! stop

end subroutine


! subroutine jacobi_transform_write(iw,jacdata)
! implicit double precision (a-h,o-z)


! type(jacobi_transform_data)       :: jacdata

! write (iw,"(I8.8)")       jacdata%n
! write (iw,"(I8.8)")       jacdata%krank
! write (iw,"(D30.20)")     jacdata%da
! write (iw,"(D30.20)")     jacdata%db
! write (iw,"(D30.20)")     jacdata%r
! write (iw,"(D30.20)")     jacdata%vals0
! write (iw,"(D30.20)")     jacdata%ts
! write (iw,"(D30.20)")     jacdata%whts
! write (iw,"(D30.20)")     jacdata%dnus
! write (iw,"(D30.20)")     jacdata%dnus0
! write (iw,"(D30.20)")     jacdata%xs
! write (iw,"(I8.8)")       jacdata%idxs
! write (iw,"(D30.20)")     jacdata%expvals

! end subroutine


! subroutine jacobi_transform_read(iw,jacdata)
! implicit double precision (a-h,o-z)
! type(jacobi_transform_data), intent(out)  :: jacdata

! read (iw,"(I8.8)")       jacdata%n
! read (iw,"(I8.8)")       jacdata%krank
! read (iw,"(D30.20)")     jacdata%da
! read (iw,"(D30.20)")     jacdata%db


! allocate(jacdata%r(jacdata%n,jacdata%krank))
! allocate(jacdata%z1(jacdata%n))
! allocate(jacdata%z2(jacdata%n))
! allocate(jacdata%z3(jacdata%n))
! allocate(jacdata%vals0(jacdata%n))
! allocate(jacdata%ts(jacdata%n))
! allocate(jacdata%whts(jacdata%n))
! allocate(jacdata%dnus(jacdata%n))
! allocate(jacdata%dnus0(jacdata%krank))
! allocate(jacdata%xs(jacdata%n))
! allocate(jacdata%idxs(jacdata%n))
! allocate(jacdata%expvals(jacdata%n,jacdata%krank))

! read (iw,"(D30.20)")     jacdata%r
! read (iw,"(D30.20)")     jacdata%vals0
! read (iw,"(D30.20)")     jacdata%ts
! read (iw,"(D30.20)")     jacdata%whts
! read (iw,"(D30.20)")     jacdata%dnus
! read (iw,"(D30.20)")     jacdata%dnus0
! read (iw,"(D30.20)")     jacdata%xs
! read (iw,"(I8.8)")       jacdata%idxs
! read (iw,"(D30.20)")     jacdata%expvals

! #ifdef USE_FFTW
! call dfftw_plan_dft_1d(jacdata%iplan,n,jacdata%z1,jacdata%z2,FFTW_BACKWARD,FFTW_ESTIMATE)
! #else
! allocate(jacdata%wsave(jacdata%n*3+15))
! call zffti(jacdata%n,jacdata%wsave)
! #endif

! end subroutine


subroutine jacobi_transform_prepare(expdata,n,jacdata)
implicit double precision (a-h,o-z)

type(jacobi_expansion_data)              :: expdata
type(jacobi_transform_data), intent(out) :: jacdata

!
!
!
!  Input parameters:
!
!  Output parameters:
!
!

double precision, allocatable :: cosvals(:,:),sinvals(:,:)
double complex                :: ima
double precision, pointer     :: arr(:)
data ima  / (0.0d0,1.0d0) /


! extract some information from the expdata structure
da      = expdata%da
db      = expdata%db
krank   = expdata%krank
nrec    = expdata%cd(1,1)


! allocate memory and setup the jacdata output structure 
jacdata%n     = n
jacdata%krank = krank
jacdata%da    = da
jacdata%db    = db
jacdata%nrec  = nrec

allocate( jacdata%ts(n),jacdata%whts(n) )
allocate( jacdata%dnus(n),jacdata%xs(n),jacdata%idxs(n),jacdata%dnus0(krank) )
allocate( jacdata%expvals(n,krank), jacdata%r(n,krank) )
allocate( jacdata%vals0(n,nrec) )
allocate( jacdata%z1(n), jacdata%z2(n), jacdata%z3(n)  )

! construct the quadrature
call jacobi_quad_mod(n,da,db,jacdata%ts,jacdata%whts)

jacdata%dnus0 = expdata%dnus0
do j=1,n
jacdata%dnus(j) = j-1
end do

jacdata%idxs = int(jacdata%ts*(n+0.0d0)/(2*pi))
jacdata%xs   = 2*pi/(n+0.0d0)*jacdata%idxs
jacdata%idxs = jacdata%idxs+1

dmax = expdata%dmax
dd = maxval(abs(jacdata%xs-jacdata%ts)*dmax)
call prin2("maximum value = ",dd)

! allocate temporary memory for interpolation
allocate( cosvals(n,krank), sinvals(n,krank) )

! perform the necessary interpolation
call jacobi_transform_intleft(expdata%chebdata1,expdata%nintsab,expdata%ab, &
  n,jacdata%ts,expdata%cosvals0,cosvals)

call jacobi_transform_intleft(expdata%chebdata1,expdata%nintsab,expdata%ab, &
  n,jacdata%ts,expdata%sinvals0,sinvals)



call jacobi_transform_intleft(expdata%chebdata2,expdata%nintscd,expdata%cd, &
  n-nrec,jacdata%dnus(nrec+1:n),transpose(expdata%r),jacdata%r(nrec+1:n,:))

jacdata%r(1:nrec,:)  = 0

jacdata%expvals      = cosvals + ima*sinvals


! scale the matrix
do j=1,krank
dnu                  = jacdata%dnus0(j)
dconst               = sqrt( (2*dnu+da+db+1)/pi )
jacdata%expvals(:,j) = dconst*jacdata%expvals(:,j) * exp(ima*(jacdata%ts-jacdata%xs)*dnu)
end do

! prepare the fft
#ifdef USE_FFTW
call dfftw_plan_dft_1d(jacdata%iplan,n,jacdata%z1,jacdata%z2,FFTW_BACKWARD,FFTW_ESTIMATE)
#else
!allocate( jacdata%wsave(3*n+15) )
!call zffti(n,jacdata%wsave)
#endif

! record the values of the polynomials of low degree
! dd                     = sqrt( (1+da+db) * gamma(1+da+db)* 1/gamma(1+da) * 1/gamma(1+db) )
! jacdata%vals0(:,1)     = dd*cos(jacdata%ts/2)**(db+0.5d0)*sin(jacdata%ts/2)**(da+0.5d0)

call  jacobi_recurrence2(jacdata%n,jacdata%ts,jacdata%nrec-1,da,db,jacdata%vals0)

deallocate(cosvals,sinvals)

end subroutine




subroutine jacobi_transform_forward(jacdata,x,y)
implicit double precision (a-h,o-z)

type(jacobi_expansion_data)            :: expdata
type(jacobi_transform_data)            :: jacdata
double precision                       :: x(jacdata%n),y(jacdata%n)

!
!
!  Input parameters:
!
!  Output parameters:
!
!
integer*8 :: iplan
double complex :: ima
data ima / (0.0d0,1.0d0) /

n          = jacdata%n
nrec       = jacdata%nrec
krank      = jacdata%krank
iplan      = jacdata%iplan
jacdata%z3 = 0

do j=1,krank
jacdata%z1 = jacdata%r(:,j)*x

#ifdef USE_FFTW
call dfftw_execute_dft(iplan, jacdata%z1, jacdata%z2)
jacdata%z3  = jacdata%z3 + jacdata%z2(jacdata%idxs)*jacdata%expvals(:,j)
#else
!call zfftb(n,jacdata%z1,jacdata%wsave)
!jacdata%z3  = jacdata%z3 + jacdata%z1(jacdata%idxs)*jacdata%expvals(:,j)
#endif
end do

! y = real(jacdata%z3)
! call dgemm('N','N',n,1,nrec,1.0d0,jacdata%vals0,n,x,n,1.0d0,y,n)
! y = y * sqrt(jacdata%whts)

y = (real(jacdata%z3) + matmul(jacdata%vals0,x(1:nrec)))*sqrt(jacdata%whts) 

end subroutine


subroutine jacobi_transform_backward(jacdata,x,y)
implicit double precision (a-h,o-z)

type(jacobi_expansion_data)            :: expdata
type(jacobi_transform_data)            :: jacdata
double precision                       :: x(jacdata%n),y(jacdata%n)
double precision, allocatable          :: x2(:)
!
!
!

double complex                         :: ima
integer*8 :: iplan

ima = (0.0d0,1.0d0)


n          = jacdata%n
krank      = jacdata%krank
nrec       = jacdata%nrec
iplan      = jacdata%iplan
jacdata%z3 = 0


ifodd = 0
if ((n/2)*2 .ne. n) ifodd = 1


! scale by the weights
allocate(x2(n))
x2 = x * sqrt(jacdata%whts)


jacdata%z2 = 0

if (ifodd .eq. 1) then

do j=1,krank
jacdata%z1        = jacdata%expvals(:,j)*x2

jacdata%z2(1:n/2) = jacdata%z1(1:n/2:2)+jacdata%z1(2:n/2:2)
jacdata%z2(n/2+1) = jacdata%z1(n)

#ifdef USE_FFTW
call dfftw_execute_dft(iplan, jacdata%z2, jacdata%z1)
jacdata%z3  = jacdata%z3 + jacdata%z1*jacdata%r(:,j)
#else
!jacdata%z1 = jacdata%z2
!call zfftb(n,jacdata%z1,jacdata%wsave)
!jacdata%z3  = jacdata%z3 + jacdata%z1*jacdata%r(:,j)
#endif

end do

else

do j=1,krank
jacdata%z1 = jacdata%expvals(:,j)*x2
jacdata%z2(1:n/2) = jacdata%z1(1:n/2:2)+jacdata%z1(2:n/2:2)

#ifdef USE_FFTW
call dfftw_execute_dft(iplan, jacdata%z2, jacdata%z1)
jacdata%z3  = jacdata%z3 + jacdata%z1*jacdata%r(:,j)
#else
!jacdata%z1 = jacdata%z2
!call zfftb(n,jacdata%z1,jacdata%wsave)
!jacdata%z3  = jacdata%z3 + jacdata%z1*jacdata%r(:,j)
#endif

end do

endif

y         = real(jacdata%z3) 
y(1:nrec) = y(1:nrec) + matmul(transpose(jacdata%vals0),x2)

end subroutine



subroutine jacobi_transform_forward_bf(da,db,k,n,ts,whts,x,y)
implicit double precision (a-h,o-z)

double precision              :: ts(n),whts(n),x(n),y(n)
!
!  If A is the (n,n) matrix representing the forward Jacobi transform, apply
!  the (k,n) submatrix A(1:k,1:n) to a user-specified vector by brute force.
!
!  Input parameters:
!    (da,db) - the parameters for the Jacobi transform
!    k - the number of rows of the submatrix of the transform to apply
!    n - the order of the transform
!    (ts,whts) - the nodes and weights of the appropriate
!       n-point Gauss-Jacobi quadrture rule
!    x - the input vector
!
!  Output parameters:
!    y - the output vector
!

double precision, allocatable :: vals(:,:)

!
!  Compute the matrix of values using the recurrence relations
!

 allocate(vals(k,n))


do i=1,k
call jacobi_recurrence(n-1,da,db,ts(i),vals(i,:))
end do

do i=1,k
vals(i,:) = vals(i,:)*sqrt(whts(i))
end do

y(1:k) = matmul(vals,x)

deallocate(vals)

end subroutine



subroutine jacobi_transform_backward_bf(da,db,k,n,ts,whts,x,y)
implicit double precision (a-h,o-z)

double precision              :: ts(n),whts(n),x(n),y(n)
!
!  If A is the (n,n) matrix representing the backward Jacobi transform, apply
!  the (k,n) submatrix A(1:k,1:n) to a user-specified vector by brute force.
!
!  Input parameters:
!    (da,db) - the parameters for the Jacobi transform
!    k - the number of rows of the submatrix of the transform to apply
!    n - the order of the transform
!    (ts,whts) - the nodes and weights of the appropriate
!       n-point Gauss-Jacobi quadrture rule
!    x - the input vector
!
!  Output parameters:
!    y - the output vector
!

double precision, allocatable :: vals(:,:)

!
!  Compute the matrix of values using the recurrence relations
!

allocate(vals(k,n))
do i=1,n
call jacobi_recurrence(k-1,da,db,ts(i),vals(:,i))
end do


y(1:k) = matmul(vals,x*sqrt(whts))

deallocate(vals)

end subroutine



subroutine jacobi_transform_intleft(chebdata,nints,ab,nts,ts,x,y)
implicit double precision (a-h,o-z)

type(chebexps_data)        :: chebdata
integer                    :: nints
double precision           :: ab(2,nints)
double precision           :: ts(nts)
double precision           :: y(:,:),x(:,:)

!
!  Left multiply a user-specified matrix by a piecewise Chebyshev interpolation matrix.
!
!    y(nts,l) = A * x(k*nints,l)
!
!
!  IMPORTANT - the list of user-specified points must be sorted in ascending order.
!
!  Input parameters:
!    chebdata - chebyshev expansion and interpolation data
!    nints - the number of intervals in the discretization scheme
!    ab - a (2,nints) array giving the endpoints of the intervals in the discretization
!      scheme
!    nts - the number of user-specified nodes
!    ts - an array containing the user specified nodes SORTED IN ASCENDING ORDER
!
!  Output parameters:
!


double precision           :: whts(chebdata%k)

eps0 = epsilon(0.0d0)
k    = chebdata%k
l    = size(x,2)
int0 = 1


do i=1,nts
t = ts(i)
do int=int0,nints-1
if (t .lt. ab(2,int)) exit
end do
a    = ab(1,int)
b    = ab(2,int)
int0 = int

xx    = (2*t - (b+a) ) /(b-a)
dsign = 1.0d0
do j=1,k
diff    = chebdata%xs(j) - xx
if (abs(diff) .lt. eps0) then
y(i,:)  = x(j+(int-1)*k,:)
goto 1000
endif
whts(j) = dsign/diff
dsign   = -dsign
end do

whts(1) = whts(1)/2
whts(k) = whts(k)/2
dd      = sum(whts)

do ll=1,l
y(i,ll)  = sum(whts*x(1+(int-1)*k:int*k,ll))/dd
end do

1000 continue

end do


end subroutine



subroutine jacobi_recurrence2(k,ts,n,da,db,vals)
implicit double precision (a-h,o-z)
integer                    :: n
double precision           :: ts(k),da,db
double precision           :: vals(k,0:n)

!
!  Evaluate the functions (1) for dnu=0,1,...,n at a specified point t
!  using the well-known recurrence relations.
!

dnu = n

do ii =1,k
t    = ts(ii)
x    = cos(t)

vals(ii,0) = sqrt(0.2d1**(-0.1d1-0.1d1*da-0.1d1*db)*(0.1d1+da+db))*sqrt(Gamma(0.1d1  &
+da+db)/(Gamma(0.1d1+da)*Gamma(0.1d1+db)))

if (n .gt. 0) then
vals(ii,1) = (sqrt(0.2d1**(-0.1d1-0.1d1*da-0.1d1*db)*(0.3d1+da+db))*(da-0.1d1*db+(2  &
+da+db)*x)*sqrt(Gamma(0.2d1+da+db)/(Gamma(0.2d1+da)*Gamma(0.2d1+  &
db))))/0.2d1
endif

do i=2,n

dd1 = (sqrt((-0.1d1+da+db+0.2d1*i)*(0.1d1+da+db+0.2d1*i))*(0.4d1*(-1+da+  &
db)*i*x+0.4d1*i**2*x+(da+db)*(da-0.1d1*db+(-2+da+  &
db)*x)))/(0.2d1*Sqrt((i*(da+i)*(db+i))/(da+db+i))*(da+db+i)*(-2+da+db+0.2d1*i))

dd2 =(2**((da+db)/0.2d1)*(-1+da+i)*(-1+db+i)*(da+db+0.2d1*i))/(i*(da+db+  &
i)*Sqrt(-3+da+db+0.2d1*i)*(-2+da+db+0.2d1*i)*Sqrt((2**(da+db)*(-1+da+  &
i)*(da+i)*(-1+db+i)*(db+i))/((-1+i)*i*(-1+da+db+i)*(da+db+i)*(1+da+db+&
0.2d1*i))))

vals(ii,i) = dd1*vals(ii,i-1) - dd2*vals(ii,i-2)


end do


!
!  Scale by r(t) now
!

rval       = 2.0d0**((1+da+db)/2) * cos(t/2)**(db+0.5d0) * sin(t/2)**(da+0.5d0)
vals(ii,:) = vals(ii,:) * rval


end do

end subroutine





subroutine jacobi_recurrence3(k,ts,n,da,db,vals)
implicit double precision (a-h,o-z)
integer                    :: n
double precision           :: ts(k),da,db
double precision           :: vals(k,0:n)

!
!  Evaluate the functions (1) for dnu=0,1,...,n at a specified point t
!  using the well-known recurrence relations.
!

dnu = n

do ii =1,k
x    = ts(ii)
! x    = cos(t)

vals(ii,0) = sqrt(0.2d1**(-0.1d1-0.1d1*da-0.1d1*db)*(0.1d1+da+db))*sqrt(Gamma(0.1d1  &
+da+db)/(Gamma(0.1d1+da)*Gamma(0.1d1+db)))

vals(ii,1) = (sqrt(0.2d1**(-0.1d1-0.1d1*da-0.1d1*db)*(0.3d1+da+db))*(da-0.1d1*db+(2  &
+da+db)*x)*sqrt(Gamma(0.2d1+da+db)/(Gamma(0.2d1+da)*Gamma(0.2d1+  &
db))))/0.2d1

do i=2,n

dd1 = (sqrt((-0.1d1+da+db+0.2d1*i)*(0.1d1+da+db+0.2d1*i))*(0.4d1*(-1+da+  &
db)*i*x+0.4d1*i**2*x+(da+db)*(da-0.1d1*db+(-2+da+  &
db)*x)))/(0.2d1*Sqrt((i*(da+i)*(db+i))/(da+db+i))*(da+db+i)*(-2+da+db+0.2d1*i))

dd2 =(2**((da+db)/0.2d1)*(-1+da+i)*(-1+db+i)*(da+db+0.2d1*i))/(i*(da+db+  &
i)*Sqrt(-3+da+db+0.2d1*i)*(-2+da+db+0.2d1*i)*Sqrt((2**(da+db)*(-1+da+  &
i)*(da+i)*(-1+db+i)*(db+i))/((-1+i)*i*(-1+da+db+i)*(da+db+i)*(1+da+db+&
0.2d1*i))))

vals(ii,i) = dd1*vals(ii,i-1) - dd2*vals(ii,i-2)

end do


!
!  Scale by r(t) now
!

t = acos(x)
rval       = 2.0d0**((1+da+db)/2) * cos(t/2)**(db+0.5d0) * sin(t/2)**(da+0.5d0)
vals(ii,:) = vals(ii,:) * rval


end do

end subroutine

end module
