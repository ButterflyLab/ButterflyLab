!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  This file contains code for constructing Gauss-Jacobi quadrature rules.  The n-point
!  Gauss-Jacobi rule associated with the parameters (da,db) is the quadrature rule
!
!         1                           n
!    \int (1-x)^da (1+x)^db p(x) ~  \sum p(x_j) w_j                                       (1)
!        -1                          j=1
! 
!  which is exact when p is a polynomial of degree 2n-1.  The n-point modified
!  Gauss-Jacoi quadrature rule is 
!
!        \pi                          
!    \int    cos(t/2)^(2db+1) sin(t/2)^(2da+1) p(cos(t)) dt ~ 
!         0                                                                               (2)
!
!                    n
!                  \sum sin(t/2)^(2db+1) cos(t/2)^(2da+1) p(cos(t_j)) w_j .
!                   j=1
!
!  It is exact when p is a polynomial of degree less than or equal to 2n-1.  Note that
!  the importance of the modified rule is that it integrates products of the functions
!  obtained by introducing the change of variables x = cos(t) into Jacobi's differential
!  equation.
!
!  The following subroutines should be regarded as publicly callable:
!
!    jacobi_quad -  construct the quadrature rule (1)
!
!    jacobi_quad_mod - construct the quadrature rule (2)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine jacobi_quad(n,da,db,xs,whts)
use chebyshev
implicit double precision (a-h,o-z)
integer, intent(in)           :: n
double precision, intent(in)  :: da,db
double precision, intent(out) :: xs(n),whts(n)

!
!  Return the n-point Gauss-Jacobi rule (1).
!
!  Input parameters:
!    n - the length of the desired quadrature rule
!    (da,db) - the parameters for Jacobi's differential equation
!
!  Output parameters:
!    xs - this user-supplied vector of length n will contain the nodes of the 
!      Gauss-Jacobi quadrature rule
!    whts - this user-supplied vector of length n will contain the quadrature
!      weights
!

double precision, allocatable :: ab(:,:),avals(:,:),psivals(:,:)
double precision, allocatable :: abinv(:,:),alphainv(:,:),alphainvp(:,:)


double precision, allocatable :: avals2(:,:),psivals2(:,:)
double precision, allocatable :: abinv2(:,:),alphainv2(:,:),alphainvp2(:,:)

type(chebexps_data)           :: chebdata
data  pi        / 3.14159265358979323846264338327950288d0 /

dnu = n
p   = dnu+(da+db+1)/2
k   = 20
call chebexps(k,chebdata)

!
!  Construct the phase function
!

nints = 50
allocate(ab(2,nints))
call jacobi_phase_disc(nints,ab)

allocate(avals(k,nints),psivals(k,nints))
allocate(avals2(k,nints),psivals2(k,nints))

call  jacobi_phase(chebdata,dnu,da,db,nints,ab,avals,psivals)
psivals = psivals*dnu
avals   = 1/avals**2


call  jacobi_phase(chebdata,dnu,db,da,nints,ab,avals2,psivals2)
psivals2 = psivals2*dnu
avals2   = 1/avals2**2


!
!  add dnu*t
!

do int=1,nints
a = ab(1,int)
b = ab(2,int)
psivals(:,int)  = psivals(:,int)  +  dnu * (chebdata%xs * (b-a)/2 + (b+a)/2)
psivals2(:,int) = psivals2(:,int) +  dnu * (chebdata%xs * (b-a)/2 + (b+a)/2)
end do


! determine the number of roots in (0,pi/2)

nroots = (psivals(k,nints/2)-psivals(1,1)) / pi

!
! Invert the phase function
!

allocate(abinv(2,nints),alphainv(k,nints),alphainvp(k,nints))

call jacobi_phase_inverse(nints,ab,chebdata%k,chebdata%xs,chebdata%aintl,chebdata%u, &
 psivals,avals,abinv,alphainv,alphainvp)


allocate(abinv2(2,nints),alphainv2(k,nints),alphainvp2(k,nints))

call jacobi_phase_inverse(nints,ab,chebdata%k,chebdata%xs,chebdata%aintl,chebdata%u, &
 psivals2,avals2,abinv2,alphainv2,alphainvp2)


!
!  Calculate the roots of P_n^(da,db)(cos(t)) in (0,pi)
!


dconst = pi *2.0d0**(da+db+1.0d0)

xs   = 0
whts = 0

xx   = pi/2 
int0 = 1
do i=1,nroots
do int=int0,nints-1
if (xx .lt. abinv(2,int)) exit
end do
a    = abinv(1,int)
b    = abinv(2,int)
int0 = int
call chebeval(a,b,k,chebdata%xs,alphainv(:,int),xx,t)
xs(n-i+1)   = t
xx          =  xx + pi
end do


xx   = pi/2 
int0 = 1
do i=1,n-nroots
do int=int0,nints-1
if (xx .lt. abinv2(2,int)) exit
end do
a    = abinv2(1,int)
b    = abinv2(2,int)
int0 = int
call chebeval(a,b,k,chebdata%xs,alphainv2(:,int),xx,t)
xs(i)   = t
xx      =  xx + pi
end do

!
!  Evaluate the corresponding weights
!


int0 = 1 
do i=1,nroots
idx = n-i+1
t = xs(idx)

do int=int0,nints-1
if (t .lt. ab(2,int)) exit
end do
a    = ab(1,int)
b    = ab(2,int)
int0 = int
call chebeval(a,b,k,chebdata%xs,avals(:,int),t,apval)
r         = cos(t/2)**(2*db+1) * sin(t/2)**(2*da+1)
whts(idx) = dconst * r/apval
 
end do



int0 = 1 
do i=1,n-nroots
idx = i
t   = xs(idx)

do int=int0,nints-1
if (t .lt. ab(2,int)) exit
end do
a    = ab(1,int)
b    = ab(2,int)
int0 = int
call chebeval(a,b,k,chebdata%xs,avals2(:,int),t,apval)
r         = cos(t/2)**(2*da+1) * sin(t/2)**(2*db+1)
whts(idx) = dconst * r/apval
xs(idx)   = pi-t 
end do

xs = cos(xs)


end subroutine



subroutine jacobi_quad_mod(n,da,db,xs,whts)
use chebyshev
implicit double precision (a-h,o-z)
integer, intent(in)           :: n
double precision, intent(in)  :: da,db
double precision, intent(out) :: xs(n),whts(n)

!
!  Return the n-point modified Gauss-Jacobi rule (2).
!
!  Input parameters:
!    n - the length of the desired quadrature rule
!    (da,db) - the parameters for Jacobi's differential equation
!
!  Output parameters:
!    xs - this user-supplied vector of length n will contain the nodes of the 
!      Gauss-Jacobi quadrature rule
!    whts - this user-supplied vector of length n will contain the quadrature
!      weights
!

double precision, allocatable :: ab(:,:),avals(:,:),psivals(:,:)
double precision, allocatable :: abinv(:,:),alphainv(:,:),alphainvp(:,:)


double precision, allocatable :: avals2(:,:),psivals2(:,:)
double precision, allocatable :: abinv2(:,:),alphainv2(:,:),alphainvp2(:,:)

type(chebexps_data)           :: chebdata
data  pi        / 3.14159265358979323846264338327950288d0 /

dnu = n
p   = dnu+(da+db+1)/2
k   = 20
call chebexps(k,chebdata)

!
!  Construct the phase function
!

nints = 50
allocate(ab(2,nints))
call jacobi_phase_disc(nints,ab)

allocate(avals(k,nints),psivals(k,nints))
allocate(avals2(k,nints),psivals2(k,nints))

call  jacobi_phase(chebdata,dnu,da,db,nints,ab,avals,psivals)
psivals = psivals*dnu
avals   = 1/avals**2


call  jacobi_phase(chebdata,dnu,db,da,nints,ab,avals2,psivals2)
psivals2 = psivals2*dnu
avals2   = 1/avals2**2


!
!  add dnu*t
!

do int=1,nints
a = ab(1,int)
b = ab(2,int)
psivals(:,int)  = psivals(:,int)  +  dnu * (chebdata%xs * (b-a)/2 + (b+a)/2)
psivals2(:,int) = psivals2(:,int) +  dnu * (chebdata%xs * (b-a)/2 + (b+a)/2)
end do


! determine the number of roots in (0,pi/2)

nroots = (psivals(k,nints/2)-psivals(1,1)) / pi

!
! Invert the phase function
!

allocate(abinv(2,nints),alphainv(k,nints),alphainvp(k,nints))

call jacobi_phase_inverse(nints,ab,chebdata%k,chebdata%xs,chebdata%aintl,chebdata%u, &
 psivals,avals,abinv,alphainv,alphainvp)


allocate(abinv2(2,nints),alphainv2(k,nints),alphainvp2(k,nints))

call jacobi_phase_inverse(nints,ab,chebdata%k,chebdata%xs,chebdata%aintl,chebdata%u, &
 psivals2,avals2,abinv2,alphainv2,alphainvp2)


!
!  Calculate the roots of P_n^(da,db)(cos(t)) in (0,pi)
!



xs   = 0
whts = 0

xx   = pi/2 
int0 = 1
do i=1,nroots
do int=int0,nints-1
if (xx .lt. abinv(2,int)) exit
end do
a    = abinv(1,int)
b    = abinv(2,int)
int0 = int
call chebeval(a,b,k,chebdata%xs,alphainv(:,int),xx,t)
xs(i)       = t
xx          =  xx + pi
end do

xx   = pi/2 
int0 = 1
do i=1,n-nroots
idx = n-i+1
do int=int0,nints-1
if (xx .lt. abinv2(2,int)) exit
end do
a    = abinv2(1,int)
b    = abinv2(2,int)
int0 = int
call chebeval(a,b,k,chebdata%xs,alphainv2(:,int),xx,t)
xs(idx)      = t
xx          =  xx + pi
end do



!
!  Evaluate the corresponding weights
!


int0 = 1 
do i=1,nroots
t = xs(i)
do int=int0,nints-1
if (t .lt. ab(2,int)) exit
end do
a    = ab(1,int)
b    = ab(2,int)
int0 = int
call chebeval(a,b,k,chebdata%xs,avals(:,int),t,apval)
!r         = cos(t/2)**(2*db+1) * sin(t/2)**(2*da+1)
whts(i)  = pi/apval
end do



int0 = 1 
do i=1,n-nroots
idx = n-i+1
t   = xs(idx)

do int=int0,nints-1
if (t .lt. ab(2,int)) exit
end do
a    = ab(1,int)
b    = ab(2,int)
int0 = int
call chebeval(a,b,k,chebdata%xs,avals2(:,int),t,apval)
whts(idx) = pi /apval
xs(idx)   = pi-t 
end do

! do i=1,n
! print *,xs(i),whts(i)
! end do
!xs = cos(xs)


end subroutine



subroutine jacobi_quad_mod2(n,da,db,xs,whts)
use chebyshev
implicit double precision (a-h,o-z)
integer, intent(in)           :: n
double precision, intent(in)  :: da,db
double precision, intent(out) :: xs(n),whts(n)

!
!  Return the modin-point Gauss-Jacobi rule (2).
!
!  Input parameters:
!    n - the length of the desired quadrature rule
!    (da,db) - the parameters for Jacobi's differential equation
!
!  Output parameters:
!    xs - this user-supplied vector of length n will contain the nodes of the 
!      Gauss-Jacobi quadrature rule
!    whts - this user-supplied vector of length n will contain the quadrature
!      weights
!

double precision, allocatable :: ab(:,:),avals(:,:),psivals(:,:)
double precision, allocatable :: abinv(:,:),alphainv(:,:),alphainvp(:,:)
type(chebexps_data)           :: chebdata
data  pi        / 3.14159265358979323846264338327950288d0 /

dnu = n
p   = dnu+(da+db+1)/2
k   = 20
call chebexps(k,chebdata)

!
!  Construct the phase function
!

nints = 50
allocate(ab(2,nints))
call jacobi_phase_disc(nints,ab)

allocate(avals(k,nints),psivals(k,nints))
call  jacobi_phase(chebdata,dnu,da,db,nints,ab,avals,psivals)
#ifdef OVERDNU
psivals = psivals*dnu
#endif
avals   = 1/avals**2

!
!  add p*t and multiply by p
!

do int=1,nints
a = ab(1,int)
b = ab(2,int)
psivals(:,int) = psivals(:,int) +  dnu * (chebdata%xs * (b-a)/2 + (b+a)/2)
end do



!
! Invert the phase function
!

allocate(abinv(2,nints),alphainv(k,nints),alphainvp(k,nints))

call jacobi_phase_inverse(nints,ab,chebdata%k,chebdata%xs,chebdata%aintl,chebdata%u, &
 psivals,avals,abinv,alphainv,alphainvp)

!
!  Calculate the roots of P_n^(da,db)(cos(t)) in (0,pi)
!

xx   = pi/2
int0 = 1
do i=1,n
do int=int0,nints-1
if (xx .lt. abinv(2,int)) exit
end do
a    = abinv(1,int)
b    = abinv(2,int)
int0 = int
call chebeval(a,b,k,chebdata%xs,alphainv(:,int),xx,xs(i))
xx = xx + pi
end do


!
!  Evaluate the corresponding weights
!

dconst = pi

int0 = 1 
do i=1,n
idx = i
t = xs(idx)

do int=int0,nints-1
if (t .lt. ab(2,int)) exit
end do
a    = ab(1,int)
b    = ab(2,int)
int0 = int
call chebeval(a,b,k,chebdata%xs,avals(:,int),t,apval)
whts(idx) = dconst/apval

end do

end subroutine




subroutine jacobi_phase_inverse(nints,ab,k,xscheb,chebintl,ucheb,alpha,alphap,&
    abinv,alphainv,alphainvp)
use chebyshev
implicit double precision (a-h,o-z)

integer, intent(in)                        :: nints,k
double precision, intent(in)               :: xscheb(k),ab(2,nints),chebintl(k,k)
double precision, intent(in)               :: alpha(k,nints),alphap(k,nints),ucheb(k,k)
double precision, intent(out)              :: alphainv(k,nints),alphainvp(k,nints),abinv(2,nints)

!
!  Compute the inverse of the phase function alpha via Newton's method.
!  
!  Input parameters:
!    (nints,ab) - the discretization scheme for representing the phase
!      function alpha
!    k - the number of terms in the piecewise Chebyshev expansions used to
!      represent the solution
!    xscheb - the nodes of the k-point Chebyshev grid on the interval [-1,1]
!    chebintl - the left Chebyshev spectral integration matrix as returned by 
!      chebexps
!    chebintr - the right Chebyshev spectral integration matrix as returned by 
!      chebexps
!    ucheb - the values-to-coefficients matrix returned by chebexps
!
!  Output parameters:
!
!    (nints,abinv) - the discretization scheme used to represent the inverse of
!      alpha
!   alphainv - a (k,nints) array specifying the values of the inverse of
!     alpha at the nodes of the k-point Chebyshev grids on the intervals
!     in the discretization scheme     
!
!   alphainvp - a (k,nints) array specifying the values of the derivative of
!     the inverse of alpha at the nodes of the k-point Chebyshev grids on the 
!     intervals in the discretization scheme        
!
!

eps0     = epsilon(0.0d0)
maxiters = 12
eps0     = sqrt(eps0)
nextra   = 3

!
!  Form the initial list of intervals for the inverse function.
!

do int=1,nints
a = ab(1,int)
b = ab(2,int)

call chebpw_eval(nints,ab,k,xscheb,alpha,a,a0)
call chebpw_eval(nints,ab,k,xscheb,alpha,b,b0)

abinv(1,int) = a0
abinv(2,int) = b0

end do

!
!  Use Newton's method to evaluate the inverse at each of the Chebyhev nodes 
!  on the grid defined by abinv; start at the right-hand side since alpha
!  is monotonically increasing.
!

do int = nints,1,-1

a0 = abinv(1,int)
b0 = abinv(2,int)
a  = ab(1,int)
b  = ab(2,int)


do i = k,1,-1

x  = (b0-a0)/2*xscheb(i) + (b0+a0)/2
t  = (b-a)/2*xscheb(i) + (b+a)/2

do iter=1,maxiters+1

if (iter .eq. maxiters+1) then
call prina("in kummer_phase_invert: maximum number of Newton iterations exceeded")
stop
endif

call chebpw_eval(nints,ab,k,xscheb,alpha,t,val)
call chebpw_eval(nints,ab,k,xscheb,alphap,t,der)

delta = (val-x)/der

if (abs(delta) .lt. eps0*(1+abs(t))) exit
t     = t - delta

end do


do iextra=1,nextra

call chebpw_eval(nints,ab,k,xscheb,alpha,t,val)
call chebpw_eval(nints,ab,k,xscheb,alphap,t,der)

delta = (val-x)/der
t     = t - delta

end do

alphainv(i,int)  = t
alphainvp(i,int) = 1.0d0/der

end do
end do

end subroutine
