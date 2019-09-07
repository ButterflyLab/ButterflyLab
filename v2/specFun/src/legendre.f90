!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  This module contains code for constructing Legendre quadratures and evaluating
!  Legendre polynomials.  
!
!  The following routines should be regarded as public:
!
!    legequad - return the nodes and weights of the n-point Legendre quadrature
!      on the interval [-1,1]
!
!    legeexps - return the nodes and weights of the n-point Legendre quadrature
!      on the interval [-1,1] as well as the matrix which maps the scaled
!      values of Legendre expansions to their coefficients
!
!    leges - evaluate the Legendre polynomials of degrees 0 through n as well 
!      as their derivatives at a specified point in the interval [-1,1]
!
!    leges0 - evaluate the Legendre polynomials of degrees 0 through n and 
!      at a specified point in the interval [-1,1]
!
!    lege - evaluate a Legendre polynomial of degree n and its derivative at  
!      a specified point on the interval [-1,1]
!
!    lege0 - evaluate a Legendre polynomial of degree n at a specified
!      point on the interval [-1,1]
!
!    legeeval - evaluate a Legendre expansion and its derivative at a 
!      specified point
!
!    legeeval0 - evaluate a Legendre expansion at a specified point
!
!    legediff - return the spectral differentiation matrix which takes
!      the scaled values of an n-term Legendre expansion at the
!      Legendre nodes to the scaled values of its derivative at the
!      Legendre nodes
!
!    legeint - return the spectral integration matrix which takes the
!      scaled values of an n-term Legendre expansion at the nodes of
!      the n-point Legnedre quadrature to the scaled values of its
!      antiderivative at those same nodes
!
!    legeq - evaluate the Legendre function of the second kind at a
!      user-specified point in the complex plane
!
!    legen - evaluate the normalized Legendre polynomials of degrees 0 through n
!      at a user-specified point in the interval [-1,1]
!
!    lege_stieltjes - evaluate P_\nu(\cos(\theta)) on the interval (1/nu, \pi/2)
!      using Stieltjes' asymptoptic formula
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module legendre

use utils

contains


subroutine legequad(n,xs,whts)
implicit double precision (a-h,o-z)
!
double precision, allocatable, intent(out) :: xs(:),whts(:)
integer n
!
!  Return the nodes and weights of the n-point Legendre quadrature on the interval
!  [-1,1].   This quadrature integrates polynomials of degree less than or equal
!  to 2n-1 exactly.
!
!  An O(n^2) algorithm is used to construct the quadrature formula when n is small
!  and the GLR algorithm is used when n is larger.
!
!  Input parameters:
!
!    n - an integer specifying the number of points to
!
!  Output parameters:
!
!    xs - upon return, a vector of length n specifying the nodes of the n-point
!      Legendre quadrature
!    whts - upon return, a of length n specifying the nodes of n-point
!      Legendre quadrature
!


double precision pols(n+1)
allocate(xs(n),whts(n))



maxiters = 10
pi       = acos(-1.0d0)

!
!  Obtain machine zero.
!

call mach_zero(eps0)

if ( n .gt. 150) then
call legequad2(n,xs,whts)
return
endif

!
!   Use Newton's method and the recurrence relation to find half of the
!   roots --- the others are obtained via symmetry.
!
!   Note that we also store the value of the derivative at each of the obtained
!   roots for use in computing the weights below.
!
ifodd = 0
nn = (n+1)/2
if (nn /= n/2) then
ifodd=1
nn=nn-1
endif

!
!  Use Tricomi's formula to  estimate the roots of P_{n+1}
!

do i =nn+1,n
   dk = i-nn-1
   dn = n
   theta = (4*(ceiling(dn/2)-dk)-1) / (4*dn+2) * pi
   x0 = 1.0d0 - (dn-1)/(8*dn**3)-1.0d0/(384.0d0*dn**4)* &
        (39.0d0-28.0d0/sin(theta)**2)
   xs(i)=x0*cos(theta)
enddo

!
!  Perform Newton iterations in order to refine the estimates.
!
do iter = 1,maxiters

!
!  Evaluate the Legendre polynomial of degree n at each point; save
!  the values of the derivative in the whts array.
!
do i=nn+1,n
call lege(n,xs(i),pols(i),whts(i))
end do

!
!  Perform one Newton iteration
!
pols(nn+1:n) = pols(nn+1:n)/whts(nn+1:n)
xs(nn+1:n)   = xs(nn+1:n) - pols(nn+1:n)

if(norm2(pols(nn+1:n)) < eps0) then
exit
endif

end do

if (iter == maxiters)  then
print *,"legequad failed!"
stop
end if

!
! Compute the weights using the derivatives we stored above.
!
do j=nn+1,n
   x       = xs(j)
   dd      = 2.0d0/(1.0d0-x**2)
   whts(j) = dd/(whts(j)**2)
end do

!
! Reflect the quadrature nodes.
!
do j=1,nn
xs(j)   = -xs(n-j+1)
whts(j) = whts(n-j+1)
end do

!
! Handle the root at 0 if n is odd.
!

if (ifodd .eq. 1) then
x0          = 0
call lege(n,x0,pol,der)
xs(nn+1)   = x0
whts(nn+1) = 2.0d0/(der**2)
endif

end subroutine


subroutine legequad2(n,xs,whts)
implicit double precision (a-h,o-z)
dimension xs(n),whts(n)
data pi /3.14159265358979323846264338327950288d0/
!
!  Use Newton's method and local Taylor expansions (i.e., the Glaser-Liu Rokhlin method)
!  to compute the n roots of the Legendre polynomial of degree n.  This routine 
!  is O(n) in the number of nodes n.

maxiters = 12
call mach_zero(eps0)
eps0=eps0*10
k        = 30
!
!  Find the roots in [0,1]
!
        ifodd = 0
        nn = (n+1)/2
        if (nn .ne. n/2) then
        ifodd=1
        nn=nn-1
        endif
!
!  Use a negative value of x0 to indicate that the procedure has
!  not been initialized.
!
!
        x0 = -1
!
        do 2000 i=nn+1,n
!
!       Use Chebyshev node as initial guess for roots ...
!
!        x1 = cos(-pi+(2*i-1)*pi/(2*n))
!
!       ... or use this somewhat better approximation.
!
        dk = i-nn-1
        dn = n
        theta = (4*(ceiling(dn/2)-dk)-1) / (4*dn+2) * pi

        x1 = 1.0d0 - (dn-1)/(8*dn**3)-1.0d0/(384.0d0*dn**4)*   &
            (39.0d0-28.0d0/sin(theta)**2)
!        
        x1=x1*cos(theta)
!
!       Conduct Newton iterations.
!
        do 2100 iter=1,maxiters
!
!       Evaluate the nth polynomial at x1 using a Taylor expansion
!       and the recurrence relation.  The first evaluation must be
!       handled specially.
!
        if (x0 .lt. 0) then
        call lege(n,x1,pol,der)
        else
        dx = x1-x0       
        call legetayl(n,k,x0,dx,pol,der)
        endif
!
!       Newton iteration.
!
        x0 = x1
        dd = -pol/der
        x1 = x0+dd
!
        if (abs(dd) .lt. eps0) then
        xs(i)=x1
        whts(i)=der
        goto 2000
        endif
!
 2100 continue
!
!       This doesn't happen.
!
        print *,"legeroots bombed!"
        stop
 2000 continue
!
!       Compute the weights using the derivatives we stored above.
!        
        do 3000 j=nn+1,n
        x       = xs(j)
        dd      = 2.0d0/(1-x**2)
        whts(j) = dd/whts(j)**2
 3000 continue
!
!       Reflect the quadrature on [-1,0] to [0,1].
!
        do 4000 j=1,nn
        xs(j)   = -xs(n-j+1)
        whts(j) = whts(n-j+1)
 4000 continue
!
!       Handle the root at 0 if n is odd.
!
        if (ifodd .eq. 1) then
        x0          = 0
        call lege(n,x0,pol,der)
!
        xs(nn+1)   = x0
        whts(nn+1) = 2.0d0/der**2
        endif
!
end subroutine


        subroutine legetayl(n,k,x,dx,pol,der)
        implicit double precision (a-h,o-z)
!
!       Evaluate the Legendre polynomial of order n using a Taylor
!       expansion and the recurrence relation for the power series
!       coefficients.
!
        k0 = min(k,n)
        p1     = pol
        p2     = der*dx
        sum    = p1+p2
        sumder = p2/dx
        do 1000 j=0,k0-2
        d = 2*x*(j+1.0d0)**2/dx*p2-(n*(n+1.0d0)-j*(j+1.0d0))*p1
        d = d/(j+1.0d0)/(j+2.0d0)*dx**2/(1.0d0-x**2)
        sum    = sum+d
        sumder = sumder+d*(j+2)/dx
        p1 = p2
        p2 = d
 1000 continue
        pol = sum
        der = sumder
        end subroutine



subroutine legeexps(n,xs,whts,u)
implicit double precision (a-h,o-z)
integer n
double precision, allocatable, intent(out) :: xs(:),whts(:),u(:,:)
!
!  Return the nodes x_1,...,x_n and weights w_1,...,w_n of the n-point
!  Legendre quadrature on the interval [-1,1].  Also, construct the
!  matrix u taking the scaled values of a n-term Legendre expansion
!  at the nodes of the n-point Legendre quadrature to its coefficients.
!  That is, if f is a polynomial of degree n-1 whose Legendre expansion is
!
!      f(x) = \sum_{j=0}^{n-1} a_j P_j(x)
!
!  then u is the matrix taking the vector
!
!     ( f(x_1) \sqrt{w_1} )
!     ( f(x_1) \sqrt{w_2} )
!             .
!             .                                                                       (1)
!             .
!     ( f(x_n) \sqrt{w_n} )
!
!   of scaled function values to the vector
!
!
!     (a_0    )
!     (a_1    )
!       .
!       .                                                                             (2)
!       .
!     (a_{n-1})
!
!
!
!
!  Input parameters:
!
!    n - the length of the Legendre quadrature to construct
!
!  Output parameters:
!
!    xs - the nodes of the n-point Legendre quadrature on the interval [-1,1]
!    whts - this array
!    u - the (n,n) matrix which maps the scaled values (1) of a
!      polynomial of degree n-1 to the coefficients (2) in its Legendre expansion
!
!      Since this matrix is orthogonal, its transpose takes the coefficients
!      of the Legendre expansion (2) to its scaled values at the nodes of
!      the n-point Legendre quadrature.
!

!
!  Call legequad to construct the n-point quadrature formula
!
call legequad(n,xs,whts)

!
!  Construct the matrix u.
!
allocate(u(n,n))

do i=1,n
call leges0(n-1,xs(i),u(:,i))
u(:,i) = u(:,i)*sqrt(whts(i))
end do

do j=1,n
dd = (2*j-1)/2.0d0
u(j,:) = u(j,:)*dd
end do

end subroutine




subroutine leges(n,x,pols,ders)
implicit double precision (a-h,o-z)

double precision, intent(out) :: pols(n+1),ders(n+1)
double precision x
integer n
!
!  Evaluate the n+1 nLegendre polynomials of degree 0 through n at the
!  point x using the standard 3-term recurrence relation.  Return the values
!  of their derivative at the point x as well.
!
!  Input parameters:
!
!    n - an integer specifying the order of polynomials which are to
!      be evaluated
!    x - the point at which the polynomials and their derivatives
!
!  Output parameters:
!
!    pols - the ith entry of this user-supplied array of length n+1
!      will contain the value of normalized Legendre polynomial of degree
!      i-1 at x
!    ders - the ith entry of this user-supplied array will contain the
!      value of the derivative of the normalized Legendre polynomial
!
!

pols(1) = 1
ders(1) = 0

if (n >= 1) then
pols(2) = x
ders(2) = 1
end if

!
!  Calculate the values of the unnormalized polynomials
!
do j=2,n
   pols(j+1)=((2.0d0*j-1.0d0)*x*pols(j)-(j-1.0d0)*pols(j-1))/j
end do

!
!  Compute the derivatives of the unnormalized polynoials
!
d=x**2.0d0-1.0d0
do j=3,n+1
   ders(j)=(j-1.0d0)*(x*pols(j)-pols(j-1))/d
end do



end subroutine


subroutine leges0(n,x,pols)
implicit double precision (a-h,o-z)

double precision, intent(out) :: pols(n+1)
double precision x
integer n
!
!  Evaluate the n+1 Legendre polynomials of degree 0 through n at the
!  point x using the standard 3-term recurrence relation.
!
!  Input parameters:
!
!    n - an integer specifying the order of polynomials which are to
!      be evaluated
!    x - the point at which the polynomials and their derivatives
!
!  Output parameters:
!
!    pols - the ith entry of this user-supplied array of length n+1
!      will contain the value of P_{i-1}(t)
!

pols(1) = 1

if (n >= 1) then
pols(2) = x
end if

!
!  Calculate the values of the polynomials
!

do j=2,n
pols(j+1)=((2.0d0*j-1.0d0)*x*pols(j)-(j-1.0d0)*pols(j-1))/j
end do


end subroutine


subroutine leges2(n,x,pols,ders,ders2)
implicit double precision (a-h,o-z)

double precision, intent(out) :: pols(n+1),ders(n+1),ders2(n+1)
double precision x
integer n
!
!  Evaluate the n+1 Legendre polynomials of degree 0 through n at the
!  point x using the standard 3-term recurrence relation.  Return the values
!  of their first and second derivatives at the point x as well.
!
!  Input parameters:
!
!    n - an integer specifying the order of polynomials which are to
!      be evaluated
!    x - the point at which the polynomials and their derivatives
!
!  Output parameters:
!
!    pols - the ith entry of this user-supplied array of length n+1
!      will contain the value of normalized Legendre polynomial of degree
!      i-1 at x
!    ders - the ith entry of this user-supplied array will contain the
!      value of the derivative of the normalized Legendre polynomial
!      of degree i-1 at x
!    ders2 - the ith entry of this user-supplied array will contain
!      the value of the second derivative of the normalized Legendre
!      polynomial of degree i-1 at x
!
!


pols(1)  = 1
ders(1)  = 0
ders2(1) = 0

if (n >= 1) then
pols(2)  = x
ders(2)  = 1
ders2(2) = 0
end if

!
!  Calculate the values of the unnormalized polynomials
!
do j=2,n
   pols(j+1)=((2.0d0*j-1.0d0)*x*pols(j)-(j-1.0d0)*pols(j-1))/j
end do

!
!  Compute the derivatives of the unnormalized polynoials
!
d=x**2.0d0-1.0d0
do j=3,n+1
   ders(j)=(j-1.0d0)*(x*pols(j)-pols(j-1))/d
end do

!
!  Compute the second derivatives of the unnormalized polynomials
!
do j=3,n+1
  ders2(j) = (j-1)/d * (pols(j)+ ( (j-1)*x-2*x) /(j-1) * ders(j) - ders(j-1))
end do



end subroutine


subroutine lege(n,x,pol,der)
implicit double precision (a-h,o-z)
integer n
double precision x,pol
!
!  Evaluate the Legendre polynomial of degree n and its derivative at
!  the point x.
!
!  Input parameters:
!
!    n - the degree of the polynomial to evaluate
!    x - the point at which the polynomial is to be evaluated
!
!  Output parameters:
!
!    pol - the value of P_n(x)
!    der - the value of P_n'(x)
!


if (n == 0) then
pol = 1
der = 0
else if (n == 1) then
pol = x
der = 1.0d0
else
p1 = 1
p2 = x

do j=2,n
   p  = ((2*j-1)*x*p2-(j-1)*p1)/j
   p1 = p2
   p2 = p
end do
!
pol = p
!
! Compute the derivative using another well-known formula.
!
der=n*(x*p2-p1)/(x**2-1.0d0)
endif

end subroutine



subroutine lege0(n,x,pol)
implicit double precision (a-h,o-z)
integer n
double precision x,pol
!
!  Evaluate the Legendre polynomial of degree n at the point x.
!
!  Input parameters:
!
!    n - the degree of the polynomial to evaluate
!    x - the point at which the polynomial is to be evaluated
!
!  Output parameters:
!
!    pol - the value of P_n(x)
!


if (n == 0) then
pol = 1.0d0
else if (n == 1) then
pol = x
else
p1 = 1
p2 = x

do j=2,n
   p  = ((2*j-1)*x*p2-(j-1)*p1)/j
   p1 = p2
   p2 = p
end do

pol = p
endif

end subroutine


subroutine legeeval(n,coefs,x,val,der)
implicit double precision (a-h,o-z)

integer n
dimension coefs(n)
double precision x,val,der

!
!  Evaluate an n-term Legendre expansion and its derivative at the point x.
!  That is, given n and a_0, a_1, ..., a_{n-1} evalute the function defined by
!
!             n
!    f(t) = \sum  a_j P_j(t),
!            j=0
!
!   where P_j is the jth Legendre polynomial, and its derivative
!   at the point x.
!
!  Input parameters:
!
!    n - the number of terms in the Legendre expansion
!    coefs - an array of length n containing the coefficients of the expansion
!    x - the point at which to evaluate the expansion
!
!  Output parameters:
!
!    val - the value of the expansion at the point x
!    der - the derivative of the expansion at the point x
!
!

dimension pols(n),ders(n)

val = 0
der = 0

call leges(n-1,x,pols,ders)


do i=1,n
val = val + pols(i)*coefs(i)
der = der + ders(i)*coefs(i)
end do


end subroutine


subroutine legeeval0(n,coefs,x,val)
implicit double precision (a-h,o-z)

integer n
dimension coefs(n)
double precision x,val

!
!  Evaluate an n-term Legendre expansion at the point x.  That is, given n 
!  and a_0, a_1, ..., a_{n-1} evalute the function defined by
!
!             n
!    f(t) = \sum  a_j P_j(t),
!            j=0
!
!   where P_j is the jth Legendre polynomial at the point x.
!
!  Input parameters:
!
!    n - the number of terms in the Legendre expansion
!    coefs - an array of length n containing the coefficients of the expansion
!    x - the point at which to evaluate the expansion
!
!  Output parameters:
!
!    val - the value of the expansion at the point x
!

dimension pols(n),ders(n)

val = 0

call leges(n-1,x,pols,ders)

do i=1,n
val = val + pols(i)*coefs(i)
end do

end subroutine


subroutine legediff(n,xs,whts,u,a)
implicit double precision (a-h,o-z)

double precision, intent(in) :: xs(n),whts(n),u(n,n)
double precision, allocatable, intent(out) :: a(:,:)
!
!  Construct the (n,n) spectral differential matrix.  That is, the matrix
!  which takes scaled values of a function f at the Legendre nodes
!  to the scaled values of its derivative at the Legendre nodes.
!
!  Input parameters:
!
!    n - the number of Legendre nodes used to represent inputs and outputs
!    (xs,whts) - the nodes and weights of the n-point Legendre quadrature
!    u - the (n,n) matrix taking the scaled values of a Legendre expansion of
!      order n-1 to its coefficients
!    v - the (n,n) matrix inverse of u
!
!  Output parameters:
!
!    a - the (n,n) spectral differentiation matrix, which is allocated
!        by the routine
!

! double precision, allocatable :: xs0(:),whts0(:)
double precision, allocatable :: b(:,:)
double precision :: pols(n),ders(n)

! call legequad(n-1,xs0,whts0)

allocate(a(n,n),b(n,n))

!
!  Evaluate the derivatives of the Legendre polynomials at the
!  Legendre nodes.
!

do i=1,n
call leges(n-1,xs(i),pols,ders)
do j=1,n
b(i,j) = ders(j) *sqrt(whts(i))
end do
end do

!
!  Multiply that matrix by the matrix which takes scaled values to
!  coefficients in order to form the output matrix.
!

a = matmul(b,u)

end subroutine



subroutine legeint(n,xs,whts,u,a)
implicit double precision (a-h,o-z)

double precision, intent(in)  :: xs(n),whts(n),u(n,n)
double precision, intent(out) :: a(n,n)
!
!  Construct the (n,n) spectral integration matrix.  That is, the matrix
!  which takes scaled values of a Legendre
!  to the scaled values of its antiderivative at Legendre nodes.
!
!  Input parameters:
!
!    n -the number of
!    (xs,whts) - the nodes and weights of the n-point Legendre quadrature
!    u - the (n,n) matrix taking the scaled values of a Legendre expansion of
!      order n-1 to its coefficients
!
!  Output parameters:
!
!    a - the (n,n) spectral integration matrix
!

double precision, allocatable :: b(:,:)
double precision :: pols(n)


allocate(b(n,n))

!
!  Use the formula
!
!
!         t
!     \int P_n(x) dx =  ( t P_n(t) - P_{n-1}(t) ) / (n+1)
!        -1
!
!
!  to form the matrix of integrals of the Legendre polynomials at the
!  Legendre nodes.
!


do i=1,n
call leges0(n-1,xs(i),pols)

b(i,1) = (1+xs(i))*sqrt(whts(i))

do j=2,n
b(i,j) = (xs(i) * pols(j) - pols(j-1))/(j)
b(i,j) = b(i,j)*sqrt(whts(i))
end do
end do


!
!  Form the spectral integration matrix.
!
a = matmul(b,u)

end subroutine


subroutine legen(n,x,pols)
implicit double precision (a-h,o-z)

double precision, intent(out) :: pols(n+1)
double precision x
integer n
!
!  Evaluate the n+1 normalized Legendre polynomials of degree 0 through n at the
!  point x using the standard 3-term recurrence relation.
!
!  Input parameters:
!
!    n - an integer specifying the order of polynomials which are to
!      be evaluated
!    x - the point at which the polynomials and their derivatives
!
!  Output parameters:
!
!    pols - the ith entry of this user-supplied array of length n+1
!      will contain the value of P_{i-1}(t)
!

pols(1) = 1

if (n >= 1) then
pols(2) = x
end if

!
!  Calculate the values of the polynomials
!
do j=2,n
   pols(j+1)=((2.0d0*j-1.0d0)*x*pols(j)-(j-1.0d0)*pols(j-1))/j
end do

!
!  Normalize the polynomials.
!

do j=1,n+1
   dd = sqrt((2*(j-1)+1)/2.0d0)
   pols(j) = pols(j)*dd
end do

end subroutine



subroutine lege_stieltjes(dnu,t,val)
implicit double precision (a-h,o-z)

double precision, intent(in)     :: dnu,t
double precision, intent(out)    :: val

data pi           /3.14159265358979323846264338327950288q0/
data sqrt2overpi / 0.797884560802865355879892119868763737q0/
data piover2     / 1.57079632679489661923132169163975144q0/

val    = 0
nterms = 17
x      = dnu+0.75d0

dd = 1-1874409467055q0/(7.0368744177664q13*x**14)+  &
7426362705q0/(1.099511627776q12*x**12)-  &
20898423q0/(8.589934592q9*x**10)+180323/(1.34217728q8*x**8)-  &
671/(5.2428799999999995q5*x**6)+21/(8.192q3*x**4)-  &
1/(6.4q1*x**2)

c0  = 1.0d0/sqrt(x) * dd
dd  = 1.0d0
dd2 = 1/sin(t)

do m=0,nterms
a   = (dnu+m+0.5d0)*t-(m+0.5d0)*piover2
val = val + c0 * cos(a)*dd
dd  = dd*dd2
c0  = c0 * (m+0.5d0)**2/(2*(m+1.0d0)*(dnu+m+1.5d0))
end do

val = val * sqrt2overpi*sqrt(dd2)

end subroutine


end module
