!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  This module contains code for representing functions via piecewise Chebyshev 
!  expansions.  More accurately, functions on an interval [a,b] can be represented either 
!  via their values  at the Chebyshev grids on a collection of subintervals [a,b] or
!  via their Chebyshev expansion coefficients on said subintervals. 
!
!  The following routines should be regarded as public:
!    
!    chebexps - construct the n-point Clenshaw-Curtis quadrature on the interval [-1,1],
!      the matrix which takes the values of a function to the Chebyshev coefficients of the 
!      polynomial interpolating the function at the nodes of the quadrature rule, as well
!      as the "left" and "right" spectral integration matrices.!
!
!    chebs - evaluate the Chebyshev polynomials of orders 0 through n at a specified
!      point in the interval [-1,1]
!
!    chebeval - evaluate a polynomial of degree n-1 whose values are given 
!      on the n-point Chebyshev grid on an interval [a,b] using the well-known 
!      barycentric interpolation formula
!
!    chebeval2 - use the Clenshaw algorithm in order to evaluate an n-term Chebyshev
!      expansion on the interval [a,b] at a specified point t in [a,b]
!
!    chebevalder - evaluate an n-term Chebyshev expansion on the interval [a,b]
!      and its derivative at a specified point x in [a,b]
!
!    chebpw_eval - evaluate a function specified by its values at k-point Chebyshev
!      grids on a collection of subintervals of [a,b] at an arbitrary point x in
!      [a,b]
!
!    chebpw_eval2 - evaluate a piecewise Chebyshev expansion given on a collection
!      of subintervals of [a,b] at a specified point x in [a,b]
!
!    chebpw_aint - apply the matrix interpolating a function from a piecewise
!      Chebyshev grid to a user-specified collection of points to a vector
!
!    chebpw_aintt - apply the transpose of the matrix interpolating a function from 
!      a piecewise Chebyshev grid to a user-specified collection of points to a vector
!
!    chebpw_aint_c - apply the matrix interpolating a function from a piecewise
!      Chebyshev grid to a user-specified collection of points to a complex-valued
!      vector
!
!    chebpw_aintt_c - apply the transpose of the matrix interpolating a function from 
!      a piecewise Chebyshev grid to a user-specified collection of points to a complex-
!      valued vector
!
!    write_chebexps - write a chebexps structure to a text file on the disk
!
!    read_chebexps - read a chebexps structure from a text file on the disk
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module chebyshev

use utils

!
!  The chebexps_data structure contains the following elements:
!
!    xs - an array containing the n quadrature nodes
!    whts - an array containing the n quadrature weights
!    u - the (n,n) matrix which takes the values of an n-term Chebyshev expansion
!      at the n quadrature nodes to the n expansion coefficients
!    v - the (n,n) matrix which takes the coefficients of an nterm-Chebyshev
!      expansion to its values at the n quadratue nodes
!    aintl - the "left" spectral integration matrix which takes the values
!      of a function f(t) on the Chebyshev nodes to the value of the function g(t) 
!      defined via the formula
!
!                     t
!          g(t) = \int  f(u) du
!                     a
!
!    aintr - the "right" spectral integration matrix which takes the values
!      of a function f(t) on the Chebyshev nodes to the value of the function g(t) 
!      defined via the formula
!
!                     t
!          g(t) = \int_ f(u) du
!                     b
!    aintr2, aintr3, aintl2, aintl3 - the higher order left and right spectral 
!      integration matrices
!

type      chebexps_data
integer                       :: k
double precision, allocatable :: xs(:),whts(:),u(:,:),v(:,:)

! left and right spectral integration matrices 
double precision, allocatable :: aintr(:,:),aintr2(:,:),aintr3(:,:)
double precision, allocatable :: aintl(:,:),aintl2(:,:),aintl3(:,:)

end type  chebexps_data


contains

subroutine chebpw_aint(chebdata,nints,ab,nts,ts,x,y)
implicit double precision (a-h,o-z)

type(chebexps_data)        :: chebdata
integer                    :: nints
double precision           :: ab(2,nints), x(chebdata%k,nints)
double precision           :: ts(nts)
double precision           :: y(nts)

!
!  Apply the matrix interpolating a function whose values are given at the
!  nodes of a piecewise Chebyshev discretization scheme to a set of user-specified 
!  points.
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
!    x - the vector to which the transpose of the matrix should be applied
!
!  Output parameters:
!    y - the output vector
!


double precision           :: whts(chebdata%k)

eps0 = epsilon(0.0d0)
int0 = 1
k    = chebdata%k

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
do l=1,k
diff    = chebdata%xs(l) - xx
if (abs(diff) .lt. eps0) then
val = x(l,int)
goto 1000
endif
whts(l) = dsign/diff
dsign   = -dsign
end do

whts(1) = whts(1)/2
whts(k) = whts(k)/2
val     = sum(whts*x(:,int))/sum(whts)

1000 continue

y(i) = val
end do


end subroutine


subroutine chebpw_aintt(chebdata,nints,ab,nts,ts,x,y)
implicit double precision (a-h,o-z)

type(chebexps_data)        :: chebdata
integer                    :: nints
double precision           :: ab(2,nints)
double precision           :: ts(nts)
double precision           :: x(nts),y(chebdata%k,nints)

!
!  Apply the transpose of the matrix interpolating a function whose values
!  are given at the nodes of a piecewise Chebyshev discretization scheme to
!  a set of user-specified points.
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
!    x - the vector to which the transpose of the matrix should be applied
!
!  Output parameters:
!    y - the output vector
!

double precision           :: whts(chebdata%k)

eps0 = epsilon(0.0d0)
int0 = 1
k    = chebdata%k
y    = 0
call prin2("x = ",x)
stop

!
!  Traverse the list of ts
!

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
do l=1,k
diff    = chebdata%xs(l) - xx
if (abs(diff) .lt. eps0) then
y(l,int) = y(l,int) + x(i)
goto 1000
endif
whts(l) = dsign/diff
dsign   = -dsign
end do
whts(1) = whts(1)/2
whts(k) = whts(k)/2

y(:,int) = y(:,int) + x(i)*whts(:)/sum(whts)
1000 continue

end do

end subroutine


subroutine chebpw_aint_c(chebdata,nints,ab,nts,ts,x,y)
implicit double precision (a-h,o-z)

type(chebexps_data)        :: chebdata
integer                    :: nints
double precision           :: ab(2,nints)
double precision           :: ts(nts)
double complex             :: y(nts), x(chebdata%k,nints)


!
!  Apply the matrix interpolating a function whose values are given at the
!  nodes of a piecewise Chebyshev discretization scheme to a set of user-specified 
!  points.
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
!    x - the vector to which the transpose of the matrix should be applied
!
!  Output parameters:
!    y - the output vector
!


double precision           :: whts(chebdata%k)

eps0 = epsilon(0.0d0)
int0 = 1
k    = chebdata%k

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
do l=1,k
diff    = chebdata%xs(l) - xx
if (abs(diff) .lt. eps0) then
val = x(l,int)
goto 1000
endif
whts(l) = dsign/diff
dsign   = -dsign
end do

whts(1) = whts(1)/2
whts(k) = whts(k)/2
val     = sum(whts*x(:,int))/sum(whts)

1000 continue

y(i) = val
end do


end subroutine


subroutine chebpw_aintt_c(chebdata,nints,ab,nts,ts,x,y)
implicit double precision (a-h,o-z)

type(chebexps_data)        :: chebdata
integer                    :: nints
double precision           :: ab(2,nints)
double precision           :: ts(nts)
double complex             :: x(nts),y(chebdata%k,nints)

!
!  Apply the transpose of the matrix interpolating a function whose values
!  are given at the nodes of a piecewise Chebyshev discretization scheme to
!  a set of user-specified points.
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
!    x - the vector to which the transpose of the matrix should be applied
!
!  Output parameters:
!    y - the output vector
!

double precision           :: whts(chebdata%k)

eps0 = epsilon(0.0d0)
int0 = 1
k    = chebdata%k
y    = 0

!
!  Traverse the list of ts
!

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
do l=1,k
diff    = chebdata%xs(l) - xx
if (abs(diff) .lt. eps0) then
y(l,int) = y(l,int) + x(i)
goto 1000
endif
whts(l) = dsign/diff
dsign   = -dsign
end do
whts(1) = whts(1)/2
whts(k) = whts(k)/2

y(:,int) = y(:,int) + x(i)*whts(:)/sum(whts)
1000 continue

end do

end subroutine





subroutine chebexps(n,chebdata)
implicit double precision (a-h,o-z)
integer                               :: n
type(chebexps_data), intent(out)      :: chebdata

!
!  Populate a chebexps structure.
!
!  Input parameters:
!    n - the number of points in the chebyshev grid
!
!  Output parameters:
!   chebdata - the chebexps_data structure containing the elements described
!    above
!   
!

double precision, allocatable :: pols(:),c(:,:),d(:,:),xx(:,:)
data pi /3.14159265358979323846264338327950288d0/

allocate(chebdata%xs(n),chebdata%whts(n),chebdata%u(n,n),chebdata%v(n,n))
allocate(pols(n+1),c(1,n),d(1,n))

chebdata%k = n
!
!  Construct the nodes
!

h = pi/(n-1)
do i=1,n
chebdata%xs(n-i+1) = cos(h*(i-1))
end do

!
!  Construct the matrix u which takes values to coefficients
!

do i=1,n
x = chebdata%xs(i)
call chebs(x,n-1,pols)
do j=1,n
chebdata%u(j,i) = pols(j)
chebdata%v(i,j) = pols(j)
end do
end do


chebdata%u(1,:) = chebdata%u(1,:)/2
chebdata%u(n,:) = chebdata%u(n,:)/2
chebdata%u(:,1) = chebdata%u(:,1)/2
chebdata%u(:,n) = chebdata%u(:,n)/2
chebdata%u = chebdata%u*2.0d0/(n-1)

!
!  Construct the weights by multiplying u^t on the left by the
!  integrals of the Chebyshev polynomials.
!

c=0
c(1,1) = 2.0d0
do i=2,n-1,2
c(1,i+1) = 1.0d0/(i+1)-1.0d0/(i-1)
end do

d = matmul(c,chebdata%u)
chebdata%whts = d(1:n,1)

!
!  Form the matrix which takes the values of a function f(t) to the values of
!
!              t
!    g(t) =  \int     f(u) du
!              a
!  

allocate(xx(n,n),chebdata%aintr(n,n),chebdata%aintl(n,n))

do i=1,n
call chebs(chebdata%xs(i),n,pols)
xx(i,1) = chebdata%xs(i)
xx(i,2) = chebdata%xs(i)**2/2.0d0
do j=3,n
xx(i,j) = 0.5d0 * (pols(j+1)/j-pols(j-1)/(j-2))
end do
end do

do i=2,n
xx(i,:) = xx(i,:) - xx(1,:)
end do
xx(1,:) = 0

chebdata%aintl = matmul(xx,chebdata%u)

!
!  Form the matrix which takes the values of a function f(t) to the values of
!
!              t
!    g(t) =  \int     f(u) du
!              b
!  

xx = 0

do i=1,n
call chebs(chebdata%xs(i),n,pols)
xx(i,1) = chebdata%xs(i)
xx(i,2) = chebdata%xs(i)**2/2.0d0
do j=3,n
xx(i,j) = 0.5d0 * (pols(j+1)/j-pols(j-1)/(j-2))
end do
end do

do i=1,n-1
xx(i,:) = xx(i,:) - xx(n,:)
end do
xx(n,:) = 0

chebdata%aintr = matmul(xx,chebdata%u)

!
!  Form the higher order spectral integration matrices
!

allocate(chebdata%aintr2(n,n),chebdata%aintr3(n,n))
chebdata%aintr2 = matmul(chebdata%aintr,chebdata%aintr)
chebdata%aintr3 = matmul(chebdata%aintr,chebdata%aintr2)

allocate(chebdata%aintl2(n,n),chebdata%aintl3(n,n))
chebdata%aintl2 = matmul(chebdata%aintl,chebdata%aintl)
chebdata%aintl3 = matmul(chebdata%aintl,chebdata%aintl2)


end subroutine



subroutine chebs(x,n,pols)
implicit double precision (a-h,o-z)

integer          :: n 
double precision :: pols(0:n),x

!
!  Evaluate the Chebyshev polynomials of degree 0 through n at a specified point
!  using the standard 3-term recurrence relation.
!
!  Input parameters:
!
!    x - point at which to evaluate the polynomials
!    n - the order of the polynomials to evaluate
!
!  Output parameters:
!
!    pols - this user-supplied and allocated array of length n+1 will
!      contain the values of the polynomials of order 0 through n upon
!      return
!

if (x .eq. 1.0d0) then
do i=0,n
pols(i) = 1.0d0
end do
return
endif

if (x .eq. -1.0d0) then
pols(0) = 1.0d0
do i=1,n
pols(i) = -pols(i-1)
end do
return
endif

pols(0) = 1.0d0
if (n .eq. 0) return

pols(1) = x
if (n .eq. 1) return

xx1 = 1.0d0
xx2 = x

do i=1,n-1
xx        = 2*x*xx2 - xx1
pols(i+1) = xx
xx1       = xx2
xx2       = xx
end do

end subroutine



subroutine chebeval(a,b,n,xs,vals,x,val)
implicit double precision (a-h,o-z)

integer :: n
double precision ::  xs(n),vals(n),x,val

!
!  Use the barycentric formula to evaluate a function given its values at the 
!  n-point Chebyshev grid on an interval [a,b].
!
!  Input parameters:
!
!    (a,b) - the interval on which the function is given
!    n - the number of nodes in the Chebyshev grid
!    xs - an array specifying the n Chevyshev node on the interval [-1,1]
!    vals - the values of the function on the n Chebyshev nodes on the
!      interval [-1,1]
!    x - the point in the interval (a,b) at which the function is to be
!      evaluated
!     
!  Output parameters:
!
!   val - the approximate value of the function at the point x
!

eps0 = epsilon(0.0d0)

xx   = (2*x - (b+a) ) /(b-a)

sum1=0
sum2=0

dd1 = 1.0d0

do i=1,n
dd=1.0d0
if (i .eq. 1 .OR. i .eq. n) dd = 0.5d0

diff = xx-xs(i)

!
!  Handle the case in which the target node coincide with one of
!  of the Chebyshev nodes.
!

if(abs(diff) .le. eps0) then
val = vals(i)
return
endif

!
!  Otherwise, construct the sums.
!

dd   = (dd1*dd)/diff
dd1  = - dd1
sum1 = sum1+dd*vals(i)
sum2 = sum2+dd
dd   = - dd
end do

val = sum1/sum2


end subroutine


subroutine chebeval2(a,b,n,coefs,x,val)
implicit double precision (a-h,o-z)

integer, intent(in)             :: n
double precision, intent(in)    :: a,b,x,coefs(n)
double precision, intent(out)   :: val

!
!  Use the Clenshaw algorithm in order to evaluate a Chebyshev expansion on the 
!  interval [a,b].
!
!  Input parameters:
!    (a,b) - the interval on which the expansion is given
!    n - the number of terms in the Chebyshev expansion
!    coefs - an array of length n specifying the expansion coefficients
!    x - the point at which to evaluate the expansion
!
!  Output parameters:
!    val - the value of the expansion at the point x
!

xx = (x - (b+a)/2.0d0) * 2.0d0/(b-a)

b2 = coefs(n)
b1 = 2*xx*b2+coefs(n-1)

do i=n-2,2,-1
b0  = coefs(i)+2*xx*b1-b2
b2 = b1
b1 = b0
end do

val = b1 * xx + (coefs(1)-b2)

end subroutine




subroutine chebpw_eval(nints,ab,k,xscheb,vals,x,val)
implicit double precision (a-h,o-z)

integer, intent(in)           :: nints,k
double precision, intent(in)  :: xscheb(k),ab(2,nints),vals(k,nints)
double precision, intent(out) :: val

!
!  Evaluate a function represented via its values at the nodes of the k-point
!  Chebyshev grids on a collection of subintervals of [a,b].
!
!  Input parameters:
!
!    (nints,ab) - arrays specifying the collection of subintervals of [a,b]
!    k - the number of terms in the Chebyshev expansions
!    xscheb - the nodes of the k-point Clenshaw-Curtis quadrature on [-1,1]
!    vals - a (k,nints) array the jth column of which gives the values of the
!      function at the nodes of the k-point Chebyshev grid in the jth
!      subinterval
!    x - the point in [a,b] at which to evaluate the function
!
!  Output parameters:
!    val - the value of the function at the point x
! 


!
!  Conduct several iterations of a binary search for the interval.
!

eps0 = epsilon(0.0d0)


niters = 7
intl   = 1
intr   = nints

do iter=1,niters
int   = (intl+intr)/2
c     = ab(1,int)
if (x .gt. c) then
intl = int
else
if (int .gt. 1) intr = int-1
endif
end do


!
!  Conduct a brute force check from here.
!

do int = intl,intr-1
b = ab(2,int)
if (x .le. b) exit
end do


!
!  Call chebeval to evaluate the expansion and the save the index
!  of the interval containing x.
!

a = ab(1,int)
b = ab(2,int)


! call chebeval(a,b,k,xscheb,vals(1,int),x,val)
! return
xx   = (2*x - (b+a) ) /(b-a)

sum1=0
sum2=0

dd1 = 1.0d0

do i=1,k
dd=1.0d0
if (i .eq. 1 .OR. i .eq. k) dd = 0.5d0

diff = xx-xscheb(i)

!
!  Handle the case in which the target node coincide with one of
!  of the Chebyshev nodes.
!

if(abs(diff) .le. eps0 ) then
val = vals(i,int)
return
endif

!
!  Otherwise, construct the sums.
!

dd   = (dd1*dd)/diff
dd1  = - dd1
sum1 = sum1+dd*vals(i,int)
sum2 = sum2+dd
dd   = - dd
end do

val = sum1/sum2


end subroutine



subroutine chebpw_eval2(nints,ab,k,coefs,x,val)
implicit double precision (a-h,o-z)

integer, intent(in)           :: nints,k
double precision, intent(in)  :: ab(2,nints),coefs(k,nints)
double precision, intent(out) :: val

!
!  Evaluate a function represented via piecewise chebyshev expansions on
!  a collection of subintervals.
!
!  Input parameters:
!    (nints,ab) - the 
!    k - an integer specifying the order of the Chebyshev expansions; on
!    
!    coefs - a (k,nints) array whose jth column specified the coefficients
!     of the function's Chebyshev expansion on the jth subinterval
!    x - the point at which to evaluate
!
!
!  Output parameters:
!    val - the value of the function at the point x
! 


double precision :: pols(k)
!
!  Conduct several iterations of a binary search for the interval.
!

niters = 6
intl   = 1
intr   = nints

do iter=1,niters
int   = (intl+intr)/2
c     = ab(1,int)
if (x .gt. c) then
intl = int
else
if (int .gt. 1) intr = int-1
endif
end do


!
!  Conduct a brute force check from here.
!

do int = intl,intr-1
b = ab(2,int)
if (x .le. b) exit
end do

a = ab(1,int)
b = ab(2,int)



!
!  Evaluate the Chebyshev expansion
!


xx = (x - (b+a)/2.0d0) * 2.0d0/(b-a)
call chebs(xx,k-1,pols)

val = 0
do i=1,k
val = val + coefs(i,int)*pols(i)
end do


return

! call chebeval2(a,b,k,coefs(1,int),x,val)


xx = (x - (b+a)/2.0d0) * 2.0d0/(b-a)
xx2 = 2*xx

b2 = coefs(k,int)
b1 = xx2*b2+coefs(k-1,int)

do i=k-2,2,-1
b0  = coefs(i,int)+xx2*b1-b2
b2 = b1
b1 = b0
end do

val = b1 * xx + (coefs(1,int)-b2)

return
end subroutine


subroutine write_chebexps(iw,chebdata)
implicit double precision (a-h,o-z)
type(chebexps_data)  :: chebdata


write (iw,"(I8.8)")       chebdata%k
write (iw,"(D30.20)")     chebdata%xs
write (iw,"(D30.20)")     chebdata%whts
write (iw,"(D30.20)")     chebdata%u
write (iw,"(D30.20)")     chebdata%v
write (iw,"(D30.20)")     chebdata%aintl
write (iw,"(D30.20)")     chebdata%aintl2
write (iw,"(D30.20)")     chebdata%aintl3
write (iw,"(D30.20)")     chebdata%aintr
write (iw,"(D30.20)")     chebdata%aintr2
write (iw,"(D30.20)")     chebdata%aintr3

end subroutine


subroutine read_chebexps(iw,chebdata)
implicit double precision (a-h,o-z)
type(chebexps_data)  :: chebdata

read (iw,"(I8.8)")       chebdata%k
k = chebdata%k
allocate(chebdata%xs(k),chebdata%whts(k))
allocate(chebdata%u(k,k),chebdata%v(k,k))
allocate(chebdata%aintl(k,k),chebdata%aintl2(k,k),chebdata%aintl3(k,k))
allocate(chebdata%aintr(k,k),chebdata%aintr2(k,k),chebdata%aintr3(k,k))
read (iw,"(D30.20)")     chebdata%xs
read (iw,"(D30.20)")     chebdata%whts
read (iw,"(D30.20)")     chebdata%u
read (iw,"(D30.20)")     chebdata%v
read (iw,"(D30.20)")     chebdata%aintl
read (iw,"(D30.20)")     chebdata%aintl2
read (iw,"(D30.20)")     chebdata%aintl3
read (iw,"(D30.20)")     chebdata%aintr
read (iw,"(D30.20)")     chebdata%aintr2
read (iw,"(D30.20)")     chebdata%aintr3

end subroutine


end module
