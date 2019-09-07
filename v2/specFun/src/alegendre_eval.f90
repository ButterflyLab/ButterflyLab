!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Copyright 2017 by James Bremer    
!
!  This program is free software: you can redistribute it and/or modify it under the terms 
!  of the GNU General Public License as published by the Free Software Foundation, either 
!  version 3 of the License, or (at your option) any later version.
!
!  This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
!  without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
!  See the GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License along with
!  this program.  If not, see <http://www.gnu.org/licenses/>.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  This file contains code for evaluating the scaled, normalized associated Legendre 
!  functions of the first and seconds kinds 
!
!
!     (  (dnu + 1/2) Gamma(dnu+dmu+1) )
!     (  ---------------------------- )^(1/2)  P_dnu^{-dmu}(\cos(t)) \sqrt{\sin(t)}       (1)
!     (        Gamma(dnu-dmu+1)       )
!
!  and 
!
!     2  (  (dnu + 1/2) Gamma(dnu+dmu+1) )
!   ---- (  ---------------------------- )^(1/2)  Q_dnu^{-dmu}(\cos(t)) \sqrt{\sin(t)}    (2)
!    Pi  (        Gamma(dnu-dmu+1)       )
!
!  on the interval (0,\pi/2).  It runs in time independent of the degree dnu and order dmu
!  and can be applied when 
!
!     0 <= dnu <= 1,000,000 and 0 <= dmu <= dnu.                                          (3)
!
!  There is also code for calculating the roots of (1) or (2) on the interval 
!  (0,\pi) which can be applied when 
!
!     10 <= dnu <= 1,000,000 and 0 <= dmu <= dnu.
!
!  The time taken to evaluate each root is independent of dnu and dmu.  
!
!  For the most part, these tasks are accomplished using precomputed tables stored in 
!  the files alegendre_data.bin1 and alegendre_data.bin2.   These tables must be read into 
!  memory by the routine  "alegendre_eval_init" before any other subroutine is called.  
!  The algorithm used for the evaluation of (1) and (2) is described in the preprint
!
!    James Bremer, "An algorithm for the numerical evaluation of the associated Legendre 
!    functions in time independent of degree and order."  arXiv:1707.03287
!
!  and the algorithm used for the calculation of roots is described in the preprint
!
!    James Bremer, "O(1) computation of the roots of the associated Legendre 
!    functions."  arXiv:???.?????
!
!  The principal subroutine for the evaluation of (1) and (2) is called  alegendre_eval.
!  We call the subset of R^3
!
!    { (dnu,dmu,t) : dnu >= 0, 0 <= dmu <= 1/2 and 0 <t <= pi }  \cup
!    { (dnu,dmu,t) : dnu >= 0, dmu > 1/2 and  tp <= t <= \pi/2 }  \cup                    (4)
!
!  where tp the turning point tp  = arcsin( sqrt(dmu^2-1/4) / (dnu+1/2) ), the
!  oscillatory region.  When (dnu,dmu,t) is in the oscillatory region, in addition 
!  to the values of (1) and (2), alegendre_eval returns the values of a nonoscillatory 
!  phase function alpha for  the associated Legendre differential equation and its 
!  derivative.  
!
!  The nonoscillatory region is 
!  
!    { (dnu,dmu,t) : dnu >=0, dmu > 1/2 and 0 < t < tp }  \cup                            (5)
!
!  When (dnu,dmu,t) is in this set, alegendre_eval returns the values of the logarithms 
!  of (1) and (2) in addition to the values of (1) and (2). 
!
!  The following routines should be regarded as publically callable:
!
!    alegendre_eval_init - this initialization routine reads the precomputed tables
!      from the files "alegendre_data.bin" and "alegendre_roots.bin" into memory.  
!      These files must be in the current working directory.  
!
!      IMPORTANT: THIS ROUTINE MUST BE CALLED BEFORE ANY OF THE OTHER SUBROUTINES
!
!    alegendre_eval - evaluate the functions (1) and (2) as well as the auxillary
!      data mentioned above
!
!    alegendre_nroots - return the number of roots of the functions (1) and (2)
!      contained in the interval (0,pi).
!      
!         ***This routine requires that dnu >= 10***
!
!    alegendre_proot - return the j^th largest root of (1) on the interval (0,\pi).
!
!         ***This routine requires that dnu >= 10***
!
!    alegendre_qroot - return the j^th largest root of (2) on the interval (0,\pi)
!
!         ***This routine requires that dnu >= 10***
!
!    alegendre_jacobi - return the j^th largest node in an ceil(n/2)-point
!      quadrature rule which integrates products of the functions
!
!    alegendre_tp - return the values of alpha, alpha' and alpha'' at the turning
!       point of the associated Legendre differential equation
!
!         ***This routine requires that dnu > 2 and dmu >= 1***
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


module alegendreeval

use besseleval

implicit double precision (a-h,o-z)

!
!  The following structure stores one of the precomputed expansions.
!
type       alegendre_expansion_data

integer                          :: ifbetas,ifover,ifsmalldmu,ifinverse
double precision                 :: dnu1,dnu2
integer                          :: ncoefs,k,nintsab,nintscd,nintsef
double precision, allocatable    :: ab(:,:),cd(:,:),ef(:,:)

!
!  Phase function and its derivative
!

integer                          :: ncoefsalpha,ncoefsalphap
double precision, allocatable    :: coefsalpha(:),coefsalphap(:)
integer, allocatable             :: iptrsalpha(:,:,:),iptrsalphap(:,:,:)

!
!  Second derivative of the phase function at the turning point
!

integer                          :: ncoefsalphapp
double precision, allocatable    :: coefsalphapp(:)
integer, allocatable             :: iptrsalphapp(:,:)

!
!  Logarithms
!

integer                          :: ncoefsbeta1,ncoefsbeta2,nintsabb
double precision, allocatable    :: coefsbeta1(:),coefsbeta2(:),abb(:,:)
integer, allocatable             :: iptrsbeta1(:,:,:),iptrsbeta2(:,:,:)

!
!  For the inverse
!

integer                          :: ncoefsalphainv,nintsinv
double precision, allocatable    :: coefsalphainv(:),abinv(:,:)
integer, allocatable             :: iptrsalphainv(:,:,:)

!
!  Some meta data regarding the precomputation
!

double precision                 :: epsrequired,epsphase,epsdisc,epscomp
double precision                 :: dmemory,dtime

!
! These are unused, except in testing
!

double precision, allocatable    :: coefsalpha0(:,:,:,:),coefsalphap0(:,:,:,:),coefsalphapp0(:,:,:)
double precision, allocatable    :: coefsbeta10(:,:,:,:),coefsbeta20(:,:,:,:)
double precision, allocatable    :: coefsalphainv0(:,:,:,:)

end type   alegendre_expansion_data



type(alegendre_expansion_data), private :: expdata1,expdata2,expdata3,expdata4
type(alegendre_expansion_data), private :: expdata5,expdata6
type(alegendre_expansion_data), private :: expdata7,expdata8
integer, private                        :: ifloaded,maxdegree
double precision, private               :: pi


data pi            / 3.14159265358979323846264338327950288d0  /
data ifloaded      / 0         /
data maxdegree     / 100000000 /

contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Main evaluation code
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine alegendre_eval(dnu,dmu,t,alpha,alphader,vallogp,vallogq,valp,valq)
implicit double precision (a-h,o-z)

!
!  Evaluate the *scaled*, *normalized* associated Legendre functions of the first 
!  and second kinds (1) and (2) of degrees nu and orders -dmu at the point cos(t), 
!  where  0 < t <= pi/2 and (dnu,dmu) is in the set (3).
!
!  When (dnu,dmu,t) is in the oscillatory region (4), also return the values of a 
!  nonoscillatory phase function alpha for the associated Legendre differential 
!  equation and its derivative.
!
!  When (dnu,dmu,t) is in the nonoscillatory region (5), also return the values of
!  the logarithms of the functions (1) and (2).
!
!  Input parameters:
!    dnu - the degree of the associated Legendre functions to evaluate
!    dmu - the order of the associated Legendre functions to evaluate
!    t - the argument at which to evaluate them
!
!  Output parameters:
!    alpha - the value of a nonoscillatory phase function if (dnu,dmu,t) is in 
!      the oscillatory region (4)
!    alphader - the value of the derivative of the nonoscillatory phase function
!      if (dnu,dmu,t) is in the oscillatory region (4)
!    vallogp - the value of the logarithm of (1) if (dnu,dmu,t) is in (5)
!    vallogq - the value of the logarithm of (2) if (dnu,dmu,t) is in (5)
!    valp - the value of (1)
!    valq - the value of (2)
!

alpha     = 0
alphader  = 0
vallogp   = 0
vallogq   = 0
valp      = 0
valq      = 0

!
!  Check to see that the precomputed expansions have been loaded
!  

if (ifloaded .eq. 0) then
print *,"alegendre_eval:  alegendre_eval_init must be called before alegendre_eval"
stop
endif

!
!  Perform range checking. 
!

if (dnu .lt. 0 .OR. dnu .gt. maxdegree   .OR. &
    dmu .lt. 0 .OR. dmu .gt. dnu         .OR. &
    t .le. 0                             .OR. &
    t .gt. pi) then
print *,"alegendre_eval: parameters out of range"
print *,"dnu = ",dnu
print *,"dmu = ",dmu
print *,"t   = ",t
stop
endif

dlambda        = dnu+0.5d0
ifoscillatory  = 1

if (dmu .gt. 0.5d0) then
deta    = sqrt(dmu**2-0.25d0)
a       = asin(deta/dlambda)
if (t .lt. a) ifoscillatory = 0
endif

ifsmalldmu = 0
if (dmu .lt. 1) ifsmalldmu = 1
call alegendre_evalabc(ifsmalldmu,dnu,dmu,a,b,c)


!
!  For small dnu, always use series expansions
!

if (dnu .lt. 2) then
call alegendre_taylor(dnu,dmu,t,alpha,alphader,vallogp,vallogq,valp,valq,ifoscillatory)
return
endif

!
!  Handle the case of 2 <= dnu < 10
!


if (dnu .lt. 10) then


if (ifsmalldmu .eq. 1) then

if (t .lt. a) then
call alegendre_taylor(dnu,dmu,t,alpha,alphader,vallogp,vallogq,valp,valq,ifoscillatory)
else
call alegendre_expeval(expdata1,dnu,dmu,t,alpha,alphader,vallogp,vallogq,valp,valq,a,b,c)
endif

else

if (t .lt. a) then
call alegendre_taylor(dnu,dmu,t,alpha,alphader,vallogp,vallogq,valp,valq,ifoscillatory)
else
call alegendre_expeval(expdata2,dnu,dmu,t,alpha,alphader,vallogp,vallogq,valp,valq,a,b,c)
endif

endif

return
endif


!
!  Now 10 <= dnu <= 10,000 
!

if (dnu .le. 10000) then

if (ifsmalldmu .eq. 1) then

if (t .lt. a) then
call alegendre_taylor(dnu,dmu,t,alpha,alphader,vallogp,vallogq,valp,valq,ifoscillatory)
else
call alegendre_expeval(expdata3,dnu,dmu,t,alpha,alphader,vallogp,vallogq,valp,valq,a,b,c)
endif

else

if (t .lt. c) then
call alegendre_taylor(dnu,dmu,t,alpha,alphader,vallogp,vallogq,valp,valq,ifoscillatory)
else
call alegendre_expeval(expdata4,dnu,dmu,t,alpha,alphader,vallogp,vallogq,valp,valq,a,b,c)
endif
endif

return
endif


!
!  Now 10,000 < dnu <= 1,000,000
!

if (dnu .le. 1000000) then


if (ifsmalldmu .eq. 1) then

if (t .lt. a) then
call alegendre_macdonald(dnu,dmu,t,vallogp,vallogq,valp,valq)
else
call alegendre_expeval(expdata5,dnu,dmu,t,alpha,alphader,vallogp,vallogq,valp,valq,a,b,c)
endif

else

if (t .lt. c) then
!call alegendre_taylor(dnu,dmu,t,alpha,alphader,vallogp,vallogq,valp,valq,ifoscillatory)
call alegendre_macdonald(dnu,dmu,t,vallogp,vallogq,valp,valq)
else
call alegendre_expeval(expdata6,dnu,dmu,t,alpha,alphader,vallogp,vallogq,valp,valq,a,b,c)
endif
endif

return
endif



!
!  Now 1,000,000 < dnu <= 100,000,000
!


if (ifsmalldmu .eq. 1) then

if (t .lt. a) then
call alegendre_macdonald(dnu,dmu,t,vallogp,vallogq,valp,valq)
else
call alegendre_expeval(expdata7,dnu,dmu,t,alpha,alphader,vallogp,vallogq,valp,valq,a,b,c)
endif

else

if (t .lt. c) then
call alegendre_macdonald(dnu,dmu,t,vallogp,vallogq,valp,valq)
else
call alegendre_expeval(expdata8,dnu,dmu,t,alpha,alphader,vallogp,vallogq,valp,valq,a,b,c)
endif
endif

return



end subroutine




subroutine alegendre_eval00(dnu,dmu,t,alphader)
implicit double precision (a-h,o-z)

!
!  Evaluate only the derivative of alpha, and only for dnu > 10 and in
!  the oscillatory regime.
!  
!


ifsmalldmu = 0
if (dmu .lt. 1) ifsmalldmu = 1


!
!  Now 10 <= dnu <= 10,000 
!

if (dnu .le. 10000) then

if (ifsmalldmu .eq. 1) then
call alegendre_expeval000(expdata3,dnu,dmu,t,alphader)
else
call alegendre_expeval000(expdata4,dnu,dmu,t,alphader)
endif

return
endif


!
!  Now 10,000 < dnu <= 1,000,000
!

if (dnu .le. 1000000) then


if (ifsmalldmu .eq. 1) then
call alegendre_expeval000(expdata5,dnu,dmu,t,alphader)
else
call alegendre_expeval000(expdata6,dnu,dmu,t,alphader)
endif

return
endif



!
!  Now 1,000,000 < dnu <= 100,000,000
!


if (ifsmalldmu .eq. 1) then

call alegendre_expeval(expdata7,dnu,dmu,t,alpha,alphader,vallogp,vallogq,valp,valq,a,b,c)
else
call alegendre_expeval(expdata8,dnu,dmu,t,alpha,alphader,vallogp,vallogq,valp,valq,a,b,c)
endif

return



end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Root evaluation code
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine alegendre_nroots(dnu,dmu,nproots,nqroots)
implicit double precision (a-h,o-z)

!
!  Return the number of roots of the functions (1) and (2) on the interval 
!  (0,\pi).
!
!  Input parameters:
!    dnu - the degree of the associated Legendre function
!    dmu - the order of the associated Legendre function
!
!  Output parameters:
!    nproots - the number of roots of (1) on the interval (0,\pi)
!    nqroots - the number of roots of (2) on the interval (0,\pi)
!
!

!
!  Compute the value of the phase function at the point pi-a, where
!  a is the left endpoint of the oscillatory region.  The value of
!  alpha at a is necessarily in the interval( -pi/2,0).
!

if (ifloaded .eq. 0) then
print *,"alegendre_root:  alegendre_eval_init must be called before this subroutine"
stop
endif

ifsmalldmu = 0
if (dmu .lt. 1) ifsmalldmu = 1

call alegendre_evalabc(ifsmalldmu,dnu,dmu,a,b,c)
call alegendre_eval(dnu,dmu,a,alpha1,alphader,vallogp,vallogq,valp,valq)
alpha2  = 2*pi + pi*(dnu-dmu)-alpha1

!
!  Roots of q occur  when \alpha = pi k while those of p occur when 
!  \alpha = pi/2 + pi k.  Count the number of such points.
!

nqroots = floor( alpha2/pi )       + 1
nproots = floor( alpha2/pi-0.5d0 ) + 1


end subroutine


subroutine alegendre_proot(dnu,dmu,j,t)
implicit double precision (a-h,o-z)

!
!  Return the j^{th} largest root of the function (1) on the interval (0,\pi).
!  For this routine, dnu >= 10.
!
!  Input parameters:
!    dnu - the degree of the associated Legendre function
!    dmu - the order of the associated Legendre function
!    j - index of the root to return
!
!  Output parameters:
!    t - the location of the specified root
!

t          = 0

if (ifloaded .eq. 0) then
print *,"alegendre_root:  alegendre_eval_init must be called before this subroutine"
stop
endif

if (dnu .lt. 10 .OR. dnu .gt. maxdegree) then
print *,"alegendre_proot: parameters out of bounds"
print *,"dnu = ",dnu,maxdegree
print *,"dmu = ",dmu
print *,"t   = ",dmu
stop
endif

x      = 2*pi + pi/2 + pi*(j-1)
xright = pi/2*(dnu-dmu) + 2*pi

!
!  Check to see if the root in in the interval (0,\pi/2) or (pi/2,pi)
!

if (x .le. xright) then
call alegendre_inverse(dnu,dmu,x,t)
else
x      = 2*xright - x
call alegendre_inverse(dnu,dmu,x,t)
t      = pi - t
endif


end subroutine


subroutine alegendre_qroot(dnu,dmu,j,t)
implicit double precision (a-h,o-z)

!
!  Return the j^{th} largest root of the function (2) on the interval (0,\pi).
!  For this routine, dnu >= 10.
!
!  Input parameters:
!    dnu - the degree of the associated Legendre function
!    dmu - the order of the associated Legendre function
!    j - index of the root to return
!
!  Output parameters:
!    t - the location of the specified root
!

t        = 0

if (ifloaded .eq. 0) then
print *,"alegendre_root:  alegendre_eval_init must be called before this subroutine"
stop
endif

if (dnu .lt. 10 .OR. dnu .gt. maxdegree) then
print *,"alegendre_qroot: parameters out of bounds"
print *,"dnu = ",dnu
print *,"dmu = ",dmu
print *,"t   = ",dmu
stop
endif

x      = 2*pi + pi*(j-1)
xright = pi*(dnu-dmu)/2+2*pi

!
!  Check to see if the root in in the interval (0,\pi/2) or (pi/2,pi)
!

if (x .le. xright) then
call alegendre_inverse(dnu,dmu,x,t)
else
x      = 2*xright-x
call alegendre_inverse(dnu,dmu,x,t)
t      = pi - t
endif

end subroutine



subroutine alegendre_jacobi(dmu,n,j,t,wht)
implicit double precision (a-h,o-z)
!
!  Return the j^th node and corresponding weight of an ceil(n/2)-point quadrature 
!  rule with nodes
!
!     t_1 < t_2 < ... < t_n
!
!  and weights
!
!     w_1, w_2, ..., w_n
!
!  such that
!
!        pi                                          n
!    \int  \sin^{2 dmu + 1}(t) p(\cos(t)) dt =  \sum  sin^{2m+1}(t) p(cos(t_j)) w_j
!         0                                         i=1
!
!  whenever p is an even polynomial of degree less than or equal to 2n-2.
!  We note that the quadrature nodes all  lie in the interval (0,\pi/2).
!
!  These quadratures are images of those described in 
!
!    Mark Tygert, "Fast algorithms for spherical harmonic expansions, III"
!    Journal of Computational Physics, 229 (2010) pp. 6181-6192.
!
!  under the change of variables x = cos(t).   They can be used in applying
!  the spherical harmonic transform.  
!
!  Input parameters:
!    dmu - the order for the quadrature rule
!    n - the number of points in the quadrature rule
!    j - the index of the node to return
!
!  Output parameters:
!    t - the required quadrature node
!    wht - the corresponding weight
!


!
!  Handle the case in which n is even
!

if ( (n/2)*2 .eq. n) then

!n0         = 2*(n/2)
n0         = n
dnu        = dmu + n0
xx         = 2*pi + pi/2 + pi*(j-1)

call alegendre_inverse(dnu,dmu,xx,t)
call alegendre_eval00(dnu,dmu,t,alphader)


dersq =  (2*dnu+1.0d0)/pi * alphader  
wht   = 2*(2*dmu+2*n0+1)/(dersq)

return
endif

!
!  Now handle the case of odd n ... when n is odd, pi/2 is one of the roots
!

n0         = 2*(n/2)+1
dnu        = dmu + n0

if (j .eq. ceiling(n/2.0d0)) then
t = pi/2
call alegendre_eval00(dnu,dmu,t,alphader)
else
xx         = 2*pi + pi/2 + pi*(j-1)
call alegendre_inverse(dnu,dmu,xx,t)
call alegendre_eval00(dnu,dmu,t,alphader)
endif

dersq =  (2*dnu+1.0d0)/pi * alphader  
wht   = 2*(2*dmu+2*n0+1)/(dersq)

end subroutine




subroutine alegendre_inverse(dnu,dmu,x,t)
implicit double precision (a-h,o-z)

!
!  Return the value of the inverse phase function at the point x.
!

if (dnu .le. 10000) then
if (dmu .lt. 1) then
call alegendre_expeval_inverse(expdata3,dnu,dmu,x,t)
else
call alegendre_expeval_inverse(expdata4,dnu,dmu,x,t)
endif
return
endif

if (dnu .le. 1000000) then
if (dmu .lt. 1) then
call alegendre_expeval_inverse(expdata5,dnu,dmu,x,t)
else
call alegendre_expeval_inverse(expdata6,dnu,dmu,x,t)
endif
return
endif


if (dmu .lt. 1) then
call alegendre_expeval_inverse(expdata7,dnu,dmu,x,t)
else
call alegendre_expeval_inverse(expdata8,dnu,dmu,x,t)
endif
return

end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Code for the evaluation of the alpha and its derivatives at the turning point
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine alegendre_tp(dnu,dmu,a,alpha,alphader,alphader2)
implicit double precision (a-h,o-z)

!
!  Return the values of the nonscillatory phase function and its first 
!  two derivatives  at the turning point. 
!
!  This routine is applicable when dnu >= 2 and 1 <= dmu <= dnu.
!
!  Input parameters:
!    dnu - degree of the associated Legendre functions 
!    dmu - order of the associated Legendre functions
!   
!  Output parameters;
!    a - the location of the turning point
!    (alpha, alphader, alphader2) - the values of alpha and its first
!      two derivatives at the point a
!

if (dmu .lt. 1 .OR. dmu .gt. dnu .OR. dnu .lt. 2 .OR. dnu .gt. 1000000) then
print *,"alegendre_tp: invalid parameter"
print *,"dnu = ",dnu
print *,"dmu = ",dmu
stop
endif

ifsmalldmu = 0
call alegendre_evalabc(ifsmalldmu,dnu,dmu,a,b,c)

if (dnu .lt. 10) then
call  alegendre_expeval00(expdata2,dnu,dmu,alpha,alphader,alphader2)
return
endif

if (dnu .lt. 10000) then
call  alegendre_expeval00(expdata4,dnu,dmu,alpha,alphader,alphader2)
return
endif

call  alegendre_expeval00(expdata6,dnu,dmu,alpha,alphader,alphader2)

end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Initialization code
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine alegendre_eval_init(dsize)
implicit double precision (a-h,o-z)

!
!  Read the precomputed expansions from the binary file "alegendre_data.bin"
!
!  Input parameters:  
!    None
!
!  Output parameters:
!    dsize - a (pretty good) estimate of the  memory occupied by the table read
!     from the disk, in megabytes
!

!
!  Do nothing and return if the table is already loaded into memory.
!

if (ifloaded .eq. 1) return

!
!  Initialize bessel_eval, which is used in alegendre_macdonald.
!

call bessel_eval_init(dsize0)

!
! Read the table into memory, set the ifloaded flag and return the
! size of the table in megabytes.
!


iw = 200
open (iw, FILE = 'alegendre_data.bin1', form = 'UNFORMATTED', status = 'OLD', &
  access = 'stream', err = 1000)
call alegendre_read_expansion(iw,expdata1)
call alegendre_read_expansion(iw,expdata2)
call alegendre_read_expansion(iw,expdata3)
call alegendre_read_expansion(iw,expdata4)
close (iw)

dsize    = expdata1%dmemory + expdata2%dmemory + expdata3%dmemory + expdata4%dmemory 

iw = 201
open (iw, FILE = 'alegendre_data.bin2', form = 'UNFORMATTED', status = 'OLD', &
  access = 'stream', err = 2000)
call alegendre_read_expansion(iw,expdata5)
call alegendre_read_expansion(iw,expdata6)
close (iw)

dsize    = dsize + expdata5%dmemory + expdata6%dmemory 

!
!  Uncomment this to allow for degrees between 1,000,000 and 100,000,000
!

iw = 202
open (iw, FILE = 'alegendre_data.bin3', form = 'UNFORMATTED', status = 'OLD', &
  access = 'stream', err = 2000)
call alegendre_read_expansion(iw,expdata7)
call alegendre_read_expansion(iw,expdata8)
close (iw)
dsize    = dsize + expdata7%dmemory + expdata8%dmemory 


ifloaded  = 1

return

1000 continue

print *,"alegendre_eval_init: unable to open and/or read alegendre_data.bin1"
stop

2000 continue

print *,"alegendre_eval_init: unable to open and/or read alegendre_data.bin2"
stop

end subroutine




subroutine alegendre_read_expansion(iw,expdata)
implicit double precision (a-h,o-z)
type(alegendre_expansion_data), intent(out) :: expdata


call alegendre_read_integer_binary(iw,expdata%ifbetas)
call alegendre_read_integer_binary(iw,expdata%ifover)
call alegendre_read_integer_binary(iw,expdata%ifsmalldmu)
call alegendre_read_integer_binary(iw,expdata%ifinverse)

call alegendre_read_double_binary(iw,expdata%dnu1)
call alegendre_read_double_binary(iw,expdata%dnu2)
call alegendre_read_double_binary(iw,expdata%dtime)
call alegendre_read_double_binary(iw,expdata%dmemory)
call alegendre_read_double_binary(iw,expdata%epsdisc)
call alegendre_read_double_binary(iw,expdata%epscomp)
call alegendre_read_double_binary(iw,expdata%epsphase)
call alegendre_read_double_binary(iw,expdata%epsrequired)


call alegendre_read_integer_binary(iw,expdata%ncoefs)
call alegendre_read_integer_binary(iw,expdata%nintsab)
call alegendre_read_integer_binary(iw,expdata%nintscd)
call alegendre_read_integer_binary(iw,expdata%nintsef)
call alegendre_read_integer_binary(iw,expdata%ncoefsalpha)
call alegendre_read_integer_binary(iw,expdata%ncoefsalphap)
call alegendre_read_integer_binary(iw,expdata%ncoefsalphapp)

allocate(expdata%ab(2,expdata%nintsab))
allocate(expdata%cd(2,expdata%nintscd))
allocate(expdata%ef(2,expdata%nintsef))
allocate(expdata%coefsalpha(expdata%ncoefsalpha))
allocate(expdata%coefsalphap(expdata%ncoefsalphap))
allocate(expdata%coefsalphapp(expdata%ncoefsalphapp))

allocate(expdata%iptrsalpha(expdata%nintsab,expdata%nintscd,expdata%nintsef))
allocate(expdata%iptrsalphap(expdata%nintsab,expdata%nintscd,expdata%nintsef))
allocate(expdata%iptrsalphapp(expdata%nintscd,expdata%nintsef))

nn = expdata%nintsab*2
call alegendre_read_double_array_binary(iw,nn,expdata%ab)
nn = expdata%nintscd*2
call alegendre_read_double_array_binary(iw,nn,expdata%cd)
nn = expdata%nintsef*2
call alegendre_read_double_array_binary(iw,nn,expdata%ef)

call alegendre_read_double_array_binary(iw,expdata%ncoefsalpha,expdata%coefsalpha)
call alegendre_read_double_array_binary(iw,expdata%ncoefsalphap,expdata%coefsalphap)
call alegendre_read_double_array_binary(iw,expdata%ncoefsalphapp,expdata%coefsalphapp)

nn = expdata%nintsab * expdata%nintscd * expdata%nintsef
call alegendre_read_integer_array_binary(iw,nn,expdata%iptrsalpha)
call alegendre_read_integer_array_binary(iw,nn,expdata%iptrsalphap)
nn = expdata%nintscd * expdata%nintsef
call alegendre_read_integer_array_binary(iw,nn,expdata%iptrsalphapp)

if (expdata%ifinverse .eq. 1) then
call alegendre_read_integer_binary(iw,expdata%nintsinv)
call alegendre_read_integer_binary(iw,expdata%ncoefsalphainv)

allocate(expdata%abinv(2,expdata%nintsinv))
allocate(expdata%coefsalphainv(expdata%ncoefsalphainv))
allocate(expdata%iptrsalphainv(expdata%nintsinv,expdata%nintscd,expdata%nintsef))

nn = expdata%nintsinv*2
call alegendre_read_double_array_binary(iw,nn,expdata%abinv)

call alegendre_read_double_array_binary(iw,expdata%ncoefsalphainv,expdata%coefsalphainv)
nn = expdata%nintsinv * expdata%nintscd * expdata%nintsef
call alegendre_read_integer_array_binary(iw,nn,expdata%iptrsalphainv)
endif

if (expdata%ifbetas .eq. 1) then
call alegendre_read_integer_binary(iw,expdata%nintsabb)
call alegendre_read_integer_binary(iw,expdata%ncoefsbeta1)
call alegendre_read_integer_binary(iw,expdata%ncoefsbeta2)

allocate(expdata%abb(2,expdata%nintsabb))
allocate(expdata%coefsbeta1(expdata%ncoefsbeta1))
allocate(expdata%coefsbeta2(expdata%ncoefsbeta2))
allocate(expdata%iptrsbeta1(expdata%nintsabb,expdata%nintscd,expdata%nintsef))
allocate(expdata%iptrsbeta2(expdata%nintsabb,expdata%nintscd,expdata%nintsef))

nn = expdata%nintsabb*2
call alegendre_read_double_array_binary(iw,nn,expdata%abb)
call alegendre_read_double_array_binary(iw,expdata%ncoefsbeta1,expdata%coefsbeta1)
call alegendre_read_double_array_binary(iw,expdata%ncoefsbeta2,expdata%coefsbeta2)
nn = expdata%nintsef *  expdata%nintscd * expdata%nintsabb
call alegendre_read_integer_array_binary(iw,nn,expdata%iptrsbeta1)
call alegendre_read_integer_array_binary(iw,nn,expdata%iptrsbeta2)
endif

end subroutine



subroutine alegendre_read_double_array_binary(iw,n,data)
implicit double precision (a-h,o-z)
double precision :: data(n)
real*8 x 
read (iw) data
end subroutine



subroutine alegendre_read_double_binary(iw,data)
implicit double precision (a-h,o-z)
real*8 :: x
read (iw) x
data = x
end subroutine


subroutine alegendre_read_double_binary16(iw,data)
implicit double precision (a-h,o-z)
real :: x
read (iw) x
data = x
end subroutine


subroutine alegendre_read_double_array_binary16(iw,n,data)
implicit double precision (a-h,o-z)
double precision :: data(n)
real :: x 

eps0 = epsilon(0.0d0)

do i=1,n
read(iw) x
data(i)=x
end do

end subroutine





subroutine alegendre_read_integer_array_binary(iw,n,idata)
implicit double precision (a-h,o-z)
integer :: idata(n)
integer*8 ix
do i=1,n
read(iw) ix
idata(i) = ix
end do
end subroutine


subroutine alegendre_read_integer_binary(iw,idata)
implicit double precision (a-h,o-z)
integer*8 i 
read (iw) i
idata = i
end subroutine



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Macdonald's asymptotic expansions for small arguments 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine alegendre_macdonald(dnu,dmu,t,vallogp,vallogq,valp,valq)
implicit double precision (a-h,o-z)

double precision :: valsj(0:15),valslogj(0:15),valsratj(0:15),dsignsj(0:15)

!
!  Evaluate the logarithms of the SCALED, NORMALIZED associated Legendre
!  functions (1) and (2) when dnu is positive and large and t is small using 
!  Macdonald's asymptotic expansions (see Proc. Royal Soc. of London, 1914 
!  pgs. 221-222).  
!
!  These expansions are valid for all 0 <= dmu <= dnu and small t.
!
!  Care has been taken to reduce the potential for numerical overflow 
!  and underflow.
!
!  Input parameters:
!    dnu - the degree of the Legendre function to evaluate; dnu >=0
!    dmu - the order of the Legendre function to evaluate; 0<= dmu <= dnu.
!    t - the argument; t must be small
!
!  Output parameters:
!    vallogp - the logarithm of the function (1)
!    vallogq - the logarithm of the function (2)
!

dimension xs(30), vals(30)

data xs  / -0.100000000000000000000000000000000000D+01,  &
           -0.959492973614497389890368057066327693D+00,  &
           -0.841253532831181168861811648919367640D+00,  &
           -0.654860733945285064056925072466293503D+00,  &
           -0.415415013001886425529274149229623161D+00,  &
           -0.142314838273285140443792668616369629D+00,  &
            0.142314838273285140443792668616369701D+00,  &
            0.415415013001886425529274149229623209D+00,  &
            0.654860733945285064056925072466293599D+00,  &
            0.841253532831181168861811648919367736D+00,  &
            0.959492973614497389890368057066327693D+00,  &
            0.100000000000000000000000000000000000D+01,  &
            0.000000000000000000000000000000000000D+00,  &
            0.000000000000000000000000000000000000D+00,  &
            0.000000000000000000000000000000000000D+00,  &
            0.000000000000000000000000000000000000D+00,  &
            0.000000000000000000000000000000000000D+00,  &
            0.000000000000000000000000000000000000D+00,  &
            0.000000000000000000000000000000000000D+00,  &
            0.000000000000000000000000000000000000D+00,  &
            0.000000000000000000000000000000000000D+00,  &
            0.000000000000000000000000000000000000D+00,  &
            0.000000000000000000000000000000000000D+00,  &
            0.000000000000000000000000000000000000D+00,  &
            0.000000000000000000000000000000000000D+00,  &
            0.000000000000000000000000000000000000D+00,  &
            0.000000000000000000000000000000000000D+00,  &
            0.000000000000000000000000000000000000D+00,  &
            0.000000000000000000000000000000000000D+00,  &
            0.000000000000000000000000000000000000D+00 /


!
!  Evaluate P
! 

call alegendre_macdonald_logp0(dnu,dmu,t,vallogp)


!
!  Now evaluate Q_dnu^{-dmu} using the connection formula and interpolation
!

diff    = 1-2*(dmu - nint(dmu))
dmu0    = nint(2*dmu)/2.0d0
diff    = abs(dmu-dmu0)


if (diff .gt. 0.005d0) then
call alegendre_macdonald_logq0(dnu,dmu,t,vallogq0,dsignq0)
!dd      = 1.0d0/cos(dmu*pi)* dsignq0  + pi/2 * tan(dmu*pi) * exp(vallogp-vallogq0)
dd      = 1.0d0/cos(dmu*pi)* dsignq0  + tan(dmu*pi) * exp(vallogp-vallogq0)
if (dd .lt. 0) dd = -dd
vallogq = vallogq0 + log(dd)
else
a = dmu0-0.1d0
b = dmu0+0.1d0
n = 12

! if (n .eq. 1) then
! xs(1) = 0.0d0
! else
! h = pi/(n-1)
! do i=1,n
! xs(n-i+1) = cos(h*(i-1))
! end do
! endif


do i=1,n
x0    = xs(i) *(b-a)/2 + (b+a)/2
call alegendre_macdonald_logq0(dnu,x0,t,vallogq0,dsignq0)
call alegendre_macdonald_logp0(dnu,x0,t,vallogp0)
!dd      = 1.0d0/cos(x0*pi)* dsignq0  + pi/2 * tan(x0*pi) * exp(vallogp0-vallogq0)
dd      = 1.0d0/cos(x0*pi)* dsignq0  + tan(x0*pi) * exp(vallogp0-vallogq0)
if (dd .lt. 0) dd = -dd

vallogq = vallogq0 + log(dd)
vals(i) = vallogq
end do

!
!  Interpolate using the barycentric formula
!

xx   = (2*dmu - (b+a) ) /(b-a)

sum1 = 0
sum2 = 0
dd1  = 1.0d0
do i=1,n
dd=1.0d0
if (i .eq. 1 .OR. i .eq. n) dd = 0.5d0
diff = xx-xs(i)

!
!  Handle the case in which the target node coincides with one of
!  of the Chebyshev nodes.
!
if(abs(diff) .le. eps0) then
val = vals(i)
exit
endif

!
!  Otherwise, calculate the sums.
!

dd   = (dd1*dd)/diff
dd1  = - dd1
sum1 = sum1+dd*vals(i)
sum2 = sum2+dd
dd   = - dd
end do

vallogq    = sum1/sum2

endif

!
!  Calculate the values of the functions
!

valp = exp(vallogp) 
valq = exp(vallogq)

end subroutine



subroutine alegendre_macdonald_logp0(dnu,dmu,t,vallogp)
implicit double precision (a-h,o-z)

!
!  Evaluate the logarithm of the NORMALIZED, SCALED associated Legendre function 
!  of the first kind (1) using Macdonald's asymptoptic expansion. 
!
!  Care has been taken to minimize the potential for numerical overflow 
!  and underflow.
!


double precision :: valsj(0:15),valslogj(0:15),valsratj(0:15),dsignsj(0:15)

x  = (2*dnu+1)   * sin(t/2)

!
!  Evaluate the ratios J_{\mu+k}(x) / J_{\mu}(x)
!


do i=0,6
call alegendre_log_besselj(dmu+i,x,valslogj(i),dsignsj(i),valsj(i))
end do


do i=0,6
valsratj(i) = dsignsj(i)*exp(valslogj(i)-valslogj(0))
end do


!
!  Compute the logarithm of P 
!

term0  = valsratj(0)


term1  = valsratj(1)/(0.2d1*x)-0.1d1*valsratj(2)+(x*valsratj(3))/0.6d1
dterm1 = (0.3d1*x*valsratj(0)-0.6d1*(1+x**2)*valsratj(1)+x*((-3+  &
x**2)*valsratj(2)+x*(0.8d1*valsratj(3)-  &
0.1d1*x*valsratj(4))))/(0.12d2*x**2)
term1  = sin(t/2)**2 * term1


term2 = (9*valsratj(2))/(0.8d1*x**2)-(29*valsratj(3))/(0.6d1*x)+  &
(31*valsratj(4))/0.12d2-(11*x*valsratj(5))/0.3d2+  &
(x**2*valsratj(6))/0.72d2
term2  = sin(t/2)**4 * term2


! term3 = (75*valsratj(3))/(0.16d2*x**3)-(751*valsratj(4))/(0.24d2*x**2)+  &
! (1381*valsratj(5))/(0.48d2*x)-(1513*valsratj(6))/0.18d3+  &
! (4943*x*valsratj(7))/0.504d4-(17*x**2*valsratj(8))/0.36d3+  &
! (x**3*valsratj(9))/0.1296d4
! term3  = sin(t/2)**6 * term3


! term4 = (3675*valsratj(4))/(0.128d3*x**4)-(20877*valsratj(5))/(0.8d2*x**3)  &
! +(98683*valsratj(6))/(0.288d3*x**2)-  &
! (110177*valsratj(7))/(0.72d3*x)+(610843*valsratj(8))/0.2016d5-  &
! (44887*x*valsratj(9))/0.1512d5+(67117*x**2*valsratj(10))/0.4536d6-  &
! (23*x**3*valsratj(11))/0.648d4+(x**4*valsratj(12))/0.31104d5
! term4  = sin(t/2)**8 * term4


! term5 = (59535*valsratj(5))/(0.256d3*x**5)-  &
! (1714703*valsratj(6))/(0.64d3*x**4)+  &
! (52827053*valsratj(7))/(0.1152d5*x**3)-  &
! (3985691*valsratj(8))/(0.144d4*x**2)+  &
! (94053343*valsratj(9))/(0.12096d6*x)-  &
! (7002137*valsratj(10))/0.6048d5+  &
! (64804643*x*valsratj(11))/0.66528d7-  &
! (643931*x**2*valsratj(12))/0.13608d7+  &
! (141703*x**3*valsratj(13))/0.108864d8-  &
! (29*x**4*valsratj(14))/0.15552d6+(x**5*valsratj(15))/0.93312d6
! term5  = sin(t/2)**10 * term5

dsum     = term0 + term1 + term2
vallogp = -dmu * log(cos(t/2)* ( dnu+0.5d0) ) + valslogj(0) + log(dsum)


!
! Apply the normalizing factor and scaling factors
!

call alegendre_gamma_ratio2(dnu,dmu,dd)
vallogp = vallogp + 0.5d0 * (log(dnu+0.5d0) + dd + log(sin(t)) )


end subroutine



subroutine alegendre_macdonald_logq0(dnu,dmu,t,vallogq,dsignq)
implicit double precision (a-h,o-z)

!
!  Evaluate the logarithm of the ABSOLUTE VALUE of the SCALED, NORMALIZED 
!  associated Legendre function 
!
!
!     2  (  (dnu + 1/2) Gamma(dnu+dmu+1) )
!   ---- (  ---------------------------- )^(1/2)  Q_dnu^{dmu}(cos(t)) \sqrt{sin(t)}       
!    Pi  (        Gamma(dnu-dmu+1)       )
!
!  when dnu is large and t is small using Macdonald's asymptotic expansion 
!  (see Proc. Royal Soc. of London, 1914 pgs. 221-222).  Also, return its sign
!
!!!!!!!!!! IMPORTANT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   THE ORDER DMU IS POSITIVE HERE!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Care has been taken to minimize the potential for numerical overflow 
!  and underflow.
!


double precision :: valsy(0:15),valslogy(0:15),valsraty(0:15),dsignsy(0:15)

x  = (2*dnu+1)   * sin(t/2)

!
!  Evaluate the ratios Y_{-\mu+k}(x) / Y_{-\mu}(x) and compute the logarithm of Q
!

do i=0,6
call alegendre_log_bessely(-dmu+i,x,valslogy(i),dsignsy(i),valsy(i))
valsraty(i) =  dsignsy(i)*exp(valslogy(i)-valslogy(0))
end do

term0 = valsraty(0)

term1 = valsraty(1)/(0.2d1*x)-0.1d1*valsraty(2)+(x*valsraty(3))/0.6d1
term1 = term1*sin(t/2)**2

term2 = (9*valsraty(2))/(0.8d1*x**2)-(29*valsraty(3))/(0.6d1*x)+  &
(31*valsraty(4))/0.12d2-(11*x*valsraty(5))/0.3d2+  &
(x**2*valsraty(6))/0.72d2
term2 = term2*sin(t/2)**4

! term3 = (75*valsraty(3))/(0.16d2*x**3)-(751*valsraty(4))/(0.24d2*x**2)+  &
! (1381*valsraty(5))/(0.48d2*x)-(1513*valsraty(6))/0.18d3+  &
! (4943*x*valsraty(7))/0.504d4-(17*x**2*valsraty(8))/0.36d3+  &
! (x**3*valsraty(9))/0.1296d4
! term3 = term3*sin(t/2)**6

! term4 =  (3675*valsraty(4))/(0.128d3*x**4)-(20877*valsraty(5))/(0.8d2*x**3)+  &
! (98683*valsraty(6))/(0.288d3*x**2)-(110177*valsraty(7))/(0.72d3*x)+  &
! (610843*valsraty(8))/0.2016d5-(44887*x*valsraty(9))/0.1512d5+  &
! (67117*x**2*valsraty(10))/0.4536d6-(23*x**3*valsraty(11))/0.648d4+  &
! (x**4*valsraty(12))/0.31104d5
! term4 = term4*sin(t/2)**8

! term5 = (59535*valsraty(5))/(0.256d3*x**5)-  &
! (1714703*valsraty(6))/(0.64d3*x**4)+  &
! (52827053*valsraty(7))/(0.1152d5*x**3)-  &
! (3985691*valsraty(8))/(0.144d4*x**2)+  &
! (94053343*valsraty(9))/(0.12096d6*x)-(7002137*valsraty(10))/0.6048d5  &
! +(64804643*x*valsraty(11))/0.66528d7-  &
! (643931*x**2*valsraty(12))/0.13608d7+  &
! (141703*x**3*valsraty(13))/0.108864d8-  &
! (29*x**4*valsraty(14))/0.15552d6+(x**5*valsraty(15))/0.93312d6
! term5 = term5 * sin(t/2)**10


dsum     = term0 + term1 + term2 
dsignq = -1.0d0
if (dsum .lt. 0) then
dsignq = -dsignq
dsum   = -dsum
endif

vallogq  = log(pi/2) + dmu * log(cos(t/2)*(dnu+0.5d0))  + valslogy(0) + log(dsum)

!
! Apply the normalizing and scaling factors
!

call alegendre_gamma_ratio2(dnu,-dmu,dd)
vallogq = vallogq + log(2.0d0/pi) + 0.5d0 * (log(dnu+0.5d0) + dd + log(sin(t)) )

end subroutine



subroutine alegendre_gamma_ratio2(dnu,dmu,val)
implicit double precision (a-h,o-z)

!
!  Evaluate log ( Gamma(dnu+dmu+1) / Gamma(dnu-dmu+1) )
!


if (dnu .lt. 10 .OR. abs(dmu) .gt. dnu/1000.0d0) then
val = log_gamma(dnu+dmu+1) - log_gamma(dnu-dmu+1)
else

x = dnu+1
y = dmu

dsum = 1-(0.1d1*y)/x-(y*(1-0.3d1*y+0.2d1*y**2))/(0.6d1*x**2)+  &
(y**2*(1-0.3d1*y+0.2d1*y**2))/(0.6d1*x**3)+(y*(6+0.5d1*y-  &
0.9d2*y**2+0.155d3*y**3-0.96d2*y**4+  &
0.2d2*y**5))/(0.36d3*x**4)-(y**2*(6+0.5d1*y-0.9d2*y**2+  &
0.155d3*y**3-0.96d2*y**4+0.2d2*y**5))/(0.36d3*x**5)-(y*(360+  &
0.126d3*y-0.2863d4*y**2-0.1323d4*y**3+0.14385d5*y**4-  &
0.18711d5*y**5+0.10518d5*y**6-0.2772d4*y**7+  &
0.28d3*y**8))/(0.4536d5*x**6)+(y**2*(360+0.126d3*y-  &
0.2863d4*y**2-0.1323d4*y**3+0.14385d5*y**4-0.18711d5*y**5+  &
0.10518d5*y**6-0.2772d4*y**7+0.28d3*y**8))/(0.4536d5*x**7)+  &
(y*(45360+0.7956d4*y-0.32274d6*y**2-0.58505d5*y**3+  &
0.8106d6*y**4+0.207158d6*y**5-0.205002d7*y**6+  &
0.2238855d7*y**7-0.115776d7*y**8+0.323336d6*y**9-  &
0.4704d5*y**10+0.28d4*y**11))/(0.54432d7*x**8)-(y**2*(45360+  &
0.7956d4*y-0.32274d6*y**2-0.58505d5*y**3+0.8106d6*y**4+  &
0.207158d6*y**5-0.205002d7*y**6+0.2238855d7*y**7-  &
0.115776d7*y**8+0.323336d6*y**9-0.4704d5*y**10+  &
0.28d4*y**11))/(0.54432d7*x**9)-(y*(5443200+0.54648d6*y-  &
0.37374084d8*y**2-0.3394248d7*y**3+0.82588957d8*y**4+  &
0.6408765d7*y**5-0.101524412d9*y**6-0.12413874d8*y**7+  &
0.154651321d9*y**8-0.150057435d9*y**9+0.72418826d8*y**10-  &
0.20401128d8*y**11+0.3409472d7*y**12-0.31416d6*y**13+  &
0.1232d5*y**14))/(0.3592512d9*x**10)



val  = 2*y*log(x) + log(dsum)

endif

end subroutine



subroutine alegendre_log_besselj(dmu,t,vallogj,dsignj,valj)
implicit double precision (a-h,o-z)
!
!  Evaluate log( | J_\dmu(t) | ) and  J_\dmu(t).
! 

if (dmu .gt. 0) then
call bessel_eval(dmu,t,alpha,alphader,vallogj,vallogy,valj,valy)
dsignj  =  1.0d0

if (vallogj .eq. 0) then
dsignj  = 1.0d0
vallogj = log(abs(valj))
if (valj .lt. 0) dsignj = -1.0d0
endif

return
endif

call bessel_eval(-dmu,t,alpha,alphader,vallogj,vallogy,valj,valy)

if (vallogj .eq. 0) then
valj   = cos(pi*dmu)*valj-sin(pi*dmu)*valy
dsignj = 1.0d0
if (valj .lt. 0) dsignj = -1.0d0
vallogj = log(abs(valj))
return
endif

dd     = cos(pi*dmu)*exp(vallogj-vallogy) - sin(pi*dmu)

if (dd .lt. 0) then
dd     = -dd
dsignj = -1.0d0
endif

vallogj = vallogyj + log ( dd ) 
valj    = dsignj*exp(vallogj)

end subroutine



subroutine alegendre_log_bessely(dmu,t,vallogy,dsigny,valy)
implicit double precision (a-h,o-z)

!
!  Evaluate log( | Y_\dmu(t) | ) and  Y_\dmu(t); dmu might be negative.
! 

if (dmu .gt. 0) then
call bessel_eval(dmu,t,alpha,alphader,vallogj,vallogy,valj,valy)
dsigny  =  -1.0d0
if (vallogy .eq. 0) then
dsigny  = 1.0d0
vallogy = log(abs(valy))
if (valy .lt. 0) dsigny = -1.0d0
endif

return
endif


call bessel_eval(-dmu,t,alpha,alphader,vallogj,vallogy,valj,valy)


dsignj0 = 1.0d0
dsigny0 = -1.0d0

if (vallogy .eq. 0) then
vallogy = log(abs(valy))
vallogj = log(abs(valj))
if (valy .gt. 0) dsigny0 = 1.0d0
if (valj .lt. 0) dsignj0 = -1.0d0
endif

dd     = -sin(pi*dmu)*dsignj0*exp(vallogj-vallogy) + cos(pi*dmu)*dsigny0

dsigny = 1.0d0
if (dd .lt. 0) then
dsigny = -1.0d0
dd     = -dd
endif

vallogy = vallogy + log ( dd ) 
valy    = dsigny*exp(vallogy)


end subroutine




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Series expansions for small arguments and parameters
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine alegendre_taylor(dnu,dmu,t,alpha,alphader,vallogp,vallogq, &
  valp,valq,ifoscillatory)
implicit double precision (a-h,o-z)

!
!  Evaluate the SCALED, NORMALIZED associated Legendre functions of the first
!  and second kinds (1) and (2) of degree dnu and order -dmu at the point cos(t)
!  using series expansions.
!

!
!  In the nonoscillatory regime, compute the logarithms.
!

if (ifoscillatory .eq. 0) then
call alegendre_ptaylor_log(dnu,dmu,t,vallogp,dsignp)
call alegendre_qtaylor_log(dnu,dmu,t,vallogq)
valp    = exp(vallogp)
valq    = exp(vallogq)
return
endif

!
!  In the oscillatory regime, use the series expansions and compute
!  alpha' and alpha as needed.  This should only be done for small
!  values of dnu, dmu and t which are less than the first zeros
!  of P and Q.
!


call alegendre_ptaylor(dnu,dmu,t,valp)
call alegendre_qtaylor(dnu,dmu,t,valq)

alphader = (2*dnu+1.0d0)/pi*1/(valp**2 + valq**2)
alpha    = atan2(-valq,valp) + 2*pi
vallogp  = 0
vallogq  = 0

end subroutine




subroutine alegendre_ptaylor(dnu,dmu,t,val)
implicit double precision (a-h,o-z)

!
!  Evaluate the SCALED, NORMALIZED associated Legendre function (1)
!  using a Taylor expansion. 
!
!  This routine should only be used in the oscillatory regime, and then
!  only for small values of dnu and dmu.
!
!  Input parameters:
!    dnu - the degree of the *normalized* associated Legende function to evaluate
!    dmu - the order of the *normalize* associated Legendre function to evaluate
!    t - the point on the interval (0,\pi/2) at which to evaluate the function
!  
!  Output parameters:
!    val - the value of P_\nu^{-\mu} (\cos(t))
!

eps0     = epsilon(0.0d0)/10.0d0
i1       = 0
i2       = 49
if (floor(dnu) .eq. dnu) then
i2 = dnu
endif

if (dmu .lt. 0 .AND. ( floor(dmu) .eq. dmu )) then
i1 = -dmu
endif


dsum = 0 
dd   = 1.0d0/alegendre_gamma(dmu+i1+1.0d0) * 1/alegendre_gamma(i1+1.0d0) * sin(t/2)**(2*i1) * (-1)**(i1)


if (i1 .ne. 0) then
dd   = dd * alegendre_gamma(dnu+i1+1.0d0)/alegendre_gamma(dnu-i1+1.0d0)
endif

xx   = sin(t/2)**2

do i=i1,i2
dsum = dsum  + dd

dd   = -dd * (dnu+i+1.0d0)*(dnu-i) * 1.0d0/(i+1.0d0) * 1.0d0/(dmu+i+1.0d0) * xx 
if (abs(dd) .lt. eps0*abs(dsum)) exit
end do

val = tan(t/2)**dmu * dsum

!
!  Now apply the normalization and scaling factors
!

call alegendre_gamma_ratio2(dnu,dmu,dd)
val = val * sqrt( (dnu+0.5d0) * exp(dd)) * sqrt(sin(t))

end subroutine




subroutine alegendre_qtaylor(dnu,dmu,t,val)
implicit double precision (a-h,o-z)

!
!  Evaluate the SCALED, NORMALIZED associated Legendre function of the
!  second kind (2) using series expansions and interpolation.
!
!  This version of the routine should only be used in the oscillatory regime
!  when dnu and dmu are both smallish.
!
!  Input parameters:
!    dnu - the degree of the *normalized* associated Legende function to evaluate
!    dmu - the order of the *normalize* associated Legendre function to evaluate
!    t - the point on the interval (0,\pi/2) at which to evaluate the function
!  
!  Output parameters:
!    val - the value of P_\nu^{-\mu} (\cos(t))
!
!  Right now, the 14-point grid is used and the points have been precomputed.

dimension xs(30), vals(30)
data xs  / -0.100000000000000000000000000000000000D+01,  &
           -0.959492973614497389890368057066327693D+00,  &
           -0.841253532831181168861811648919367640D+00,  &
           -0.654860733945285064056925072466293503D+00,  &
           -0.415415013001886425529274149229623161D+00,  &
           -0.142314838273285140443792668616369629D+00,  &
            0.142314838273285140443792668616369701D+00,  &
            0.415415013001886425529274149229623209D+00,  &
            0.654860733945285064056925072466293599D+00,  &
            0.841253532831181168861811648919367736D+00,  &
            0.959492973614497389890368057066327693D+00,  &
            0.100000000000000000000000000000000000D+01,  &
            0.000000000000000000000000000000000000D+00,  &
            0.000000000000000000000000000000000000D+00,  &
            0.000000000000000000000000000000000000D+00,  &
            0.000000000000000000000000000000000000D+00,  &
            0.000000000000000000000000000000000000D+00,  &
            0.000000000000000000000000000000000000D+00,  &
            0.000000000000000000000000000000000000D+00,  &
            0.000000000000000000000000000000000000D+00,  &
            0.000000000000000000000000000000000000D+00,  &
            0.000000000000000000000000000000000000D+00,  &
            0.000000000000000000000000000000000000D+00,  &
            0.000000000000000000000000000000000000D+00,  &
            0.000000000000000000000000000000000000D+00,  &
            0.000000000000000000000000000000000000D+00,  &
            0.000000000000000000000000000000000000D+00,  &
            0.000000000000000000000000000000000000D+00,  &
            0.000000000000000000000000000000000000D+00,  &
            0.000000000000000000000000000000000000D+00 /


eps0  = epsilon(0.0d0)
mu    = nint(dmu)
diff  = abs(dmu-mu)


!
!  When dmu is not close to an integer, use the well-known formula (which
!  is simplified somewhat by normalization)
!

if (diff .gt. 0.01d0) then
call alegendre_qtaylor0(dnu,dmu,t,val)
return
endif


!
!  Otherwise, interpolate in the dmu variable
!

a = mu-0.1d0
b = mu+0.1d0
n = 12


!
!  The commented code constructs the n-point Chebyshev grid on the interval [-1,1].
!

! if (n .ne. 12) then
! if (n .eq. 1) then
! xs(1) = 0.0d0
! else
! h = pi/(n-1)
! do i=1,n
! xs(n-i+1) = cos(h*(i-1))
! end do
! endif
! endif

!
!  Evaluate Y_\nu at the nodes of the Chebyshev grid on the interval [a,b]
!

do i=1,n
x0    = xs(i) *(b-a)/2 + (b+a)/2
call alegendre_qtaylor0(dnu,x0,t,vals(i))
end do

!
!  Interpolate using the barycentric formula
!

xx   = (2*dmu - (b+a) ) /(b-a)
sum1 = 0
sum2 = 0
dd1  = 1.0d0
do i=1,n
dd=1.0d0
if (i .eq. 1 .OR. i .eq. n) dd = 0.5d0
diff = xx-xs(i)

!
!  Handle the case in which the target node coincides with one of
!  of the Chebyshev nodes.
!
if(abs(diff) .le. eps0) then
val = vals(i)
return
endif

!
!  Otherwise, calculate the sums.
!

dd   = (dd1*dd)/diff
dd1  = - dd1
sum1 = sum1+dd*vals(i)
sum2 = sum2+dd
dd   = - dd
end do

val    = sum1/sum2

end subroutine


subroutine alegendre_qtaylor0(dnu,dmu,t,val)
implicit double precision (a-h,o-z)

!
!  Evaluate the SCALED, NORMALIZED associated Legendre function (1)
!  using Taylor expansions in the event that mu is not close to an integer.
!

call alegendre_ptaylor(dnu,dmu,t,val1)
call alegendre_ptaylor(dnu,-dmu,t,val2)

!val = 1.0d0/(sin(dmu*pi)) * (val2 - val1*cos(dmu*pi))

val = (val2/sin(pi*dmu) - val1*cos(dmu*pi)/sin(dmu*pi))

end subroutine




subroutine alegendre_ptaylor_log(dnu,dmu,t,vallog,dsign)
implicit double precision (a-h,o-z)

!
!  Calculate the value of the NORMALIZED, SCALED, associated legendre
!  function of the first kind (1).
!
!  Input parameters:
!    dnu - the degree of the *normalized* associated Legende function to evaluate
!    dmu - the order of the *normalize* associated Legendre function to evaluate
!    t - the point on the interval (0,\pi/2) at which to evaluate the function
!  
!  Output parameters:
!    dsign - the sign of P_nu^{-mu}(cos(t))
!    vallog - logarithm of the absolute value of P_dnu^{-dmu} ( cos (t ))
!

eps0     = epsilon(0.0d0)/10.0d0
i1       = 0
i2       = 49

if (floor(dnu) .eq. dnu) then
i2 = dnu
endif

if (dmu .lt. 0 .AND. ( floor(dmu) .eq. dmu )) then
i1 = -dmu
endif

dsign = (-1.0d0)**i1
dd   = -log_gamma(dmu+i1+1.0d0)-log_gamma(i1+1.0d0) + 2*i1 * log(sin(t/2)) 


if (i1 .ne. 0) then
dd   = dd + log_gamma(dnu+i1+1.0d0) - log_gamma(dnu-i1+1.0d0)
endif


ddd = alegendre_gamma(dmu+i1+1.0d0)

if (ddd .lt. 0)  dsign = -dsign

dd2   = 1.0d0
dsum  = 0 
xx    = sin(t/2)**2


do i=i1,i2
dsum = dsum  + dd2
dd2   = -dd2 * (dnu+i+1.0d0)*(dnu-i) * 1.0d0/(i+1.0d0) * 1.0d0/(dmu+i+1.0d0) * xx 
if (abs(dd2) .le. eps0*abs(dsum)) exit
end do

if (dsum .lt. 0) then
dsum  = -dsum
dsign  = -dsign
endif


vallog = dd + log(dsum) + dmu * log(tan(t/2))


!
!  Now the normalization factor
!

call alegendre_gamma_ratio2(dnu,dmu,dd)


vallog = vallog + 0.5d0 * (log(dnu+0.5d0) + dd + log(sin(t)) )
val    = dsign*exp(vallog)


end subroutine




subroutine alegendre_qtaylor_log(dnu,dmu,t,vallog)
implicit double precision (a-h,o-z)

!
!  Evaluate the SCALED, NORMALIZED associated Legendre function of the
!  second kind (2) using series expansions.
!
!  Input parameters:
!    dnu - the degree of the *normalized* associated Legende function to evaluate
!    dmu - the order of the *normalize* associated Legendre function to evaluate
!    t - the point on the interval (0,\pi/2) at which to evaluate the function
!  
!  Output parameters:
!    val - the value of P_\nu^{-\mu} (\cos(t))
!
!  Right now, the 14-point grid is used and the points have been precomputed.

dimension xs(30), vals(30)
data xs  / -0.100000000000000000000000000000000000D+01,  &
           -0.959492973614497389890368057066327693D+00,  &
           -0.841253532831181168861811648919367640D+00,  &
           -0.654860733945285064056925072466293503D+00,  &
           -0.415415013001886425529274149229623161D+00,  &
           -0.142314838273285140443792668616369629D+00,  &
            0.142314838273285140443792668616369701D+00,  &
            0.415415013001886425529274149229623209D+00,  &
            0.654860733945285064056925072466293599D+00,  &
            0.841253532831181168861811648919367736D+00,  &
            0.959492973614497389890368057066327693D+00,  &
            0.100000000000000000000000000000000000D+01,  &
            0.000000000000000000000000000000000000D+00,  &
            0.000000000000000000000000000000000000D+00,  &
            0.000000000000000000000000000000000000D+00,  &
            0.000000000000000000000000000000000000D+00,  &
            0.000000000000000000000000000000000000D+00,  &
            0.000000000000000000000000000000000000D+00,  &
            0.000000000000000000000000000000000000D+00,  &
            0.000000000000000000000000000000000000D+00,  &
            0.000000000000000000000000000000000000D+00,  &
            0.000000000000000000000000000000000000D+00,  &
            0.000000000000000000000000000000000000D+00,  &
            0.000000000000000000000000000000000000D+00,  &
            0.000000000000000000000000000000000000D+00,  &
            0.000000000000000000000000000000000000D+00,  &
            0.000000000000000000000000000000000000D+00,  &
            0.000000000000000000000000000000000000D+00,  &
            0.000000000000000000000000000000000000D+00,  &
            0.000000000000000000000000000000000000D+00 /


eps0  = epsilon(0.0d0)
mu    = nint(dmu)
diff  = abs(dmu-mu)


!
!  When dmu is not close to an integer, use the well-known formula (which
!  is simplified somewhat by normalization)
!



if (diff .gt. 0.005d0) then
call alegendre_qtaylor0_log(dnu,dmu,t,vallog,dsign)
return
else

!
!  Otherwise, interpolate in the dmu variable
!

a = mu-0.1d0
b = mu+0.1d0
n = 12


!
!  The commented code constructs the n-point Chebyshev grid on the interval [-1,1].
!  Right now, the 14-point grid is used and the points have been precomputed.
!

! if (n .eq. 1) then
! xs(1) = 0.0d0
! else
! h = pi/(n-1)
! do i=1,n
! xs(n-i+1) = cos(h*(i-1))
! end do
! endif

!
!  Evaluate Y_\nu at the nodes of the Chebyshev grid on the interval [a,b]
!

do i=1,n
x0    = xs(i) *(b-a)/2 + (b+a)/2
call alegendre_qtaylor0_log(dnu,x0,t,vals(i),dsign)
end do

!
!  Interpolate using the barycentric formula
!

xx   = (2*dmu - (b+a) ) /(b-a)
sum1 = 0
sum2 = 0
dd1  = 1.0d0
do i=1,n
dd=1.0d0
if (i .eq. 1 .OR. i .eq. n) dd = 0.5d0
diff = xx-xs(i)

!
!  Handle the case in which the target node coincides with one of
!  of the Chebyshev nodes.
!
if(abs(diff) .le. eps0) then
val = vals(i)
return
endif

!
!  Otherwise, calculate the sums.
!

dd   = (dd1*dd)/diff
dd1  = - dd1
sum1 = sum1+dd*vals(i)
sum2 = sum2+dd
dd   = - dd
end do

vallog  = sum1/sum2
return

endif


end subroutine




subroutine alegendre_qtaylor0_log(dnu,dmu,t,vallog,dsign)
implicit double precision (a-h,o-z)

!
!  Calculate vallog and dsign such that |dsign|=1 and
!  
!    Q_nu^{-\mu}(cos(t)) = dsign * exp ( vallog )
!  
!

call alegendre_ptaylor_log(dnu,dmu,t,vallog1,dsign1)
call alegendre_ptaylor_log(dnu,-dmu,t,vallog2,dsign2)

dd    = 1.0d0/(sin(dmu*pi)) * (1 - dsign2*dsign1*exp(vallog1-vallog2)*cos(dmu*pi))


dsign = 1.0d0
if (dd .lt. 0) then
dd   = -dd
dsign = 1.0d0
endif

vallog = vallog2 + log(dd)

end subroutine



function alegendre_gamma(x)
implicit double precision (a-h,o-z)
double precision :: alegendre_gamma

!
!  This routine exists to deal with a bug in certain versions of the gfortran
!  library.  The bug causes the sign of the gamma function of negative values 
!  to be calculated incorrectly in certain cases.
!

dgamma = gamma(x)
if (x .lt. 0) then
nn1 = floor(x)
if ((nn1/2)*2 .eq. nn1 .AND. dgamma .lt. 0) then
dgamma = -dgamma
endif
endif

alegendre_gamma = dgamma

end function



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Code for using a series expansion to evaluate alpha' near the right endpoint pi/2
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine alegendre_nearright(dnu,dmu,t,aval,apval,valp,valq)
implicit double precision (a-h,o-z)

!
!  Evaluate alpha and alpha' near the right endpoit pi/2 using a series expansion
!  for alpha.
!  
!  Input parameters:
!    dnu - the degree of associated Legendre function to evaluate
!    dmu - the order of the associated Legendre function to evaluate
!    t - the point at which to evaluate alpha and alpha' 
!
!  Output parameters:
!    aval - the value of alpha at the point t
!    apval - the value of alpha' at the point t
!

double precision :: as(0:15)

pi    = acos(-1.0d0)

dlambdasq = (dnu+0.5d0)**2
detasq    = dmu**2-0.25d0
delta     = t-pi/2

call alegendre_alphap0(dnu,dmu,apval)

as(0) = pi/2*(dnu-dmu)
as(1) = apval
as(2) = 0.0d0
as(3) = 0.2q1*(-detasq+dlambdasq)*as(1)-0.2q1*as(1)**3

as(5)=0.2q1*dlambdasq*as(3)-0.6q1*as(1)**2*as(3)+(0.3q1*as(3)**2)/as(1)-  &
0.2q1*detasq*(2*as(1)+as(3))

as(7)=(-0.4q1*detasq*as(1)**5*(16*as(1)+12*as(3)+as(5))+  &
0.2q1*as(1)**3*(-36*as(1)**3*as(3)**2-18*as(3)**3+2*dlambdasq*as(1)**2*as(5)-  &
6*as(1)**4*as(5)+12*as(1)*as(3)*as(5)))/(0.2q1*as(1)**5)

as(9)=-(0.4q1*detasq*as(1)**7*(272*as(1)+240*as(3)+30*as(5)+as(7))+  &
0.2q1*as(1)**4*(-270*as(3)**4+180*as(1)**4*as(3)*as(5)+225*as(1)*as(3)**2*as(5)+  &
6*as(1)**5*as(7)+2*as(1)**3*(90*as(3)**3-dlambdasq*as(7))-  &
3*as(1)**2*(10*as(5)**2+6*as(3)*as(7))))/(0.2q1*as(1)**7)

as(11)=-(0.4q1*detasq*as(1)**9*(7936*as(1)+7616*as(3)+1120*as(5)+56*as(7)+as(9))+  &
0.2q1*as(1)**5*(7560*as(3)**5-7560*as(1)*as(3)**3*as(5)+  &
84*as(1)*as(3)**2*(30*as(1)**3*as(5)+7*as(1)*as(7))+  &
12*as(1)**2*as(3)*(140*as(5)**2+28*as(1)**3*as(7)-2*as(1)*as(9))+  &
as(1)**3*(420*as(1)**2*as(5)**2-168*as(5)*as(7)-2*dlambdasq*as(1)*as(9)+  &
6*as(1)**3*as(9))))/(0.2q1*as(1)**9)

as(13)=-(0.4q1*detasq*as(1)**11*(353792*as(1)+357120*as(3)+57120*as(5)+3360*as(7)+  &
90*as(9)+as(11))+0.2q1*as(1)**6*(-340200*as(3)**6+396900*as(1)*as(3)**4*as(5)-  &
30240*as(1)**2*as(3)**3*as(7)+135*as(1)**2*as(3)**2*(-910*as(5)**2+  &
56*as(1)**3*as(7)+9*as(1)*as(9))+15*as(1)**3*as(3)*(1260*as(1)**2*as(5)**2+  &
924*as(5)*as(7)+36*as(1)**3*as(9)-2*as(1)*as(11))+2*as(1)**3*(3150*as(5)**3+  &
180*as(5)*(7*as(1)**3*as(7)-as(1)*as(9))+as(1)*(-189*as(7)**2-  &
dlambdasq*as(1)*as(11)+3*as(1)**3*as(11)))))/(0.2q1*as(1)**11)

as(15)=-(0.4q1*detasq*as(1)**13*(22368256*as(1)+23350272*as(3)+3928320*as(5)+  &
251328*as(7)+7920*as(9)+132*as(11)+as(13))+0.2q1*as(1)**7*(22453200*as(3)**7-  &
29937600*as(1)*as(3)**5*as(5)+2245320*as(1)**2*as(3)**4*as(7)-  &
89100*as(1)**2*as(3)**3*(-133*as(5)**2+as(1)*as(9))+  &
198*as(1)**3*as(3)**2*(-7140*as(5)*as(7)+90*as(1)**3*as(9)+11*as(1)*as(11))+  &
as(1)**4*(110880*as(5)**2*as(7)+396*as(1)**3*(14*as(7)**2+15*as(5)*as(9))-  &
33*as(1)*(72*as(7)*as(9)+20*as(5)*as(11))+6*as(1)**4*as(13)+  &
as(1)**2*(69300*as(5)**3-2*dlambdasq*as(13)))+36*as(1)**3*as(3)*(-34650*as(5)**3  &
+1155*as(1)*as(5)*(4*as(1)**2*as(7)+as(9))+as(1)*(1155*as(7)**2+  &
22*as(1)**3*as(11)-as(1)*as(13)))))/(0.2q1*as(1)**13)

aval   = as(0)
dd     = delta
do i=0,7
aval   = aval + as(2*i+1) * dd
dd     = dd*delta*delta*1.0d0/(2*i+3.0d0)*1.0d0/(2*i+2.0d0)
end do

apval  = 0
dd     = 1.0d0
do i=0,7
apval   = apval + as(2*i+1) * dd
dd     = dd*delta*delta*1.0d0/(2*i+2.0d0)*1.0d0/(2*i+1.0d0)
end do


dd   =  sqrt((1+2*dnu)/(pi*apval))
valp =  cos(aval ) * dd
valq = -sin(aval ) * dd

end subroutine


subroutine alegendre_alphap0(dnu,dmu,val)
implicit double precision (a-h,o-z)

!
!  This routine calculates the value of the derivative of the nonoscillatory
!  phase function for the associated Legendre equation (1) at the point 0  
!  using an asymptotic expansions and the obvious recurrence relations.
!
!  The value of \alpha'(0) is equal to
!
!     - (2 dnu + 1)                        1
!     --------------  --------------------------------------------
!           pi                (P_nu^mu(0)^2 + (2/pi Q_nu^mu(0))^2
!
!  which is equal to (see, for instance, HTF Volume 1, p. 145)
!
!      2 Gamma ( 1/2 (dnu - dmu + 2 ) ) * Gamma ( 1/2 * (dnu + dmu + 2) )
!      ------------------------------------------------------------------
!        Gamma ( 1/2 (dnu - dmu + 1 ) ) * Gamma ( 1/2 * (dnu + dmu + 1) )
!
!  This routine runs in time independent of dmu and dnu and achieves roughly 
!  30 digit accracy when run in extended precision.
!
!  Input parameters:
!    dnu - degree 
!    dmu - order
!
!  Output parameters:
!    val - the value of the nonoscillatory phase function for the associated 
!      Legendre equation (1) at the point 0
!

x = (dnu-dmu+1)/2
y = (dnu+dmu+1)/2

call alegendre_gammratio1(x,val1)
call alegendre_gammratio1(y,val2)

val  = 2*val1*val2

end subroutine


subroutine alegendre_gammratio1(x,val)
implicit double precision (a-h,o-z)
!
!  Evaluate the ratio Gamma(x+1/2) / Gamma(x) either directly (for small x)
!  or using an asymptoptic expansion (for large x)
!

!
!  Handle the case of small values of x
!

if (x .lt. 10) then
val = alegendre_gamma(x+0.5d0)/alegendre_gamma(x)
return
endif

!
!  Use an asymptotic expansion for large values of x.  The function
!
!    g(y) = \Gamma(x + 1/4 +1/2) / Gamma(x + 1/4)
!
!  has a nicre expansion than Gamma(x+1/2)/Gamma(x), so we use that.
!

y = x-0.25d0

val = 0.1q1-0.838312462819535446586854867773q1/y**20-  &
0.120547316354150874426610235964q0/y**19+  &
0.966554052957437020538105450131q0/y**18+  &
0.174548668221481564590729165798q-1/y**17-  &
0.140046033565018861557528728756q0/y**16-  &
0.326950390425616299605593439992q-2/y**15+  &
0.26258284398827602679205028835q-1/y**14+  &
0.822517646662745960384199861437q-3/y**13-  &
0.661636761914508042536908760667q-2/y**12-  &
0.292163815629464806988835334778q-3/y**11+  &
0.235649039314012043178081512451q-2/y**10+  &
0.156118825543671846389770507812q-3/y**9-  &
0.12691752053797245025634765625q-2/y**8-  &
0.1556575298309326171875q-3/y**7+  &
0.12209415435791015625q-2/y**6-0.26702880859375q-4/y**5-  &
0.885009765625q-3/y**4-0.68359375q-2/y**3+0.390625q-1/y**2-  &
0.125q0/y

val = val * sqrt(x)

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Evaluation routines for the precomputed expansions
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine alegendre_expeval(expdata,dnu,dmu,t,aval,apval,bval1,bval2,valp,valq,a,b,c)
implicit double precision (a-h,o-z)

double precision, intent(in)       :: dnu,t
type(alegendre_expansion_data)     :: expdata
double precision, intent(out)      :: aval,apval,bval1,bval2

aval   = 0
apval  = 0
bval   = 0


k           = expdata%k
ifsmalldmu  = expdata%ifsmalldmu
ifover      = expdata%ifover
ifbetas     = expdata%ifbetas


if (ifover .eq. 1) then
dnu0   = 1/dnu
else
dnu0   = dnu
endif



!call alegendre_evalabc(ifsmalldmu,dnu,dmu,a,b,c)
call compute_dmu0(ifsmalldmu,dnu,dmu,dmu0)



if (t .ge. a) then

u = (t-a)/(b-a)
call alegendre_findint(expdata%nintsef,expdata%ef,dnu0,intef,e0,f0)
call alegendre_findint(expdata%nintscd,expdata%cd,dmu0,intcd,c0,d0)
call alegendre_findint(expdata%nintsab,expdata%ab,u,intab,a0,b0)


iptr1  = expdata%iptrsalpha(intab,intcd,intef)
iptr2  = expdata%iptrsalphap(intab,intcd,intef)

call alegendre_tensor_eval2(expdata%ncoefsalpha,expdata%coefsalpha,expdata%ncoefsalphap, &
  expdata%coefsalphap,iptr1,iptr2,a0,b0,c0,d0,e0,f0,u,dmu0,dnu0,aval,apval)


if (ifover .eq. 1) then
aval  = aval  * dnu
apval = apval * dnu
endif

dd   =  sqrt((1+2*dnu)/(pi*apval))
valp =  cos(aval ) * dd
valq = -sin(aval ) * dd


else

u = (t-c)/(a-c)


call alegendre_findint(expdata%nintsef,expdata%ef,dnu0,intef,e0,f0)
call alegendre_findint(expdata%nintscd,expdata%cd,dmu0,intcd,c0,d0)
call alegendre_findint(expdata%nintsabb,expdata%abb,u,intabb,a0,b0)


iptr1  = expdata%iptrsbeta1(intabb,intcd,intef)
iptr2  = expdata%iptrsbeta2(intabb,intcd,intef)

call alegendre_tensor_eval2(expdata%ncoefsbeta1,expdata%coefsbeta1, &
  expdata%ncoefsbeta2,expdata%coefsbeta2,iptr1,iptr2, &
  a0,b0,c0,d0,e0,f0,u,dmu0,dnu0,bval1,bval2)

if (ifover .eq. 1) then
bval1  = bval1  * dnu
bval2  = bval2  * dnu
endif

bval1 = bval1+dnu
bval2 = bval2-dnu

valp = exp(bval1)
valq = exp(bval2)

endif

end subroutine



subroutine alegendre_expeval000(expdata,dnu,dmu,t,apval)
implicit double precision (a-h,o-z)

type(alegendre_expansion_data) :: expdata

k           = expdata%k
ifsmalldmu  = expdata%ifsmalldmu
ifover      = expdata%ifover


if (ifover .eq. 1) then
dnu0   = 1/dnu
else
dnu0   = dnu
endif



call alegendre_evalabc(ifsmalldmu,dnu,dmu,a,b,c)
call compute_dmu0(ifsmalldmu,dnu,dmu,dmu0)



u = (t-a)/(b-a)
call alegendre_findint(expdata%nintsef,expdata%ef,dnu0,intef,e0,f0)
call alegendre_findint(expdata%nintscd,expdata%cd,dmu0,intcd,c0,d0)
call alegendre_findint(expdata%nintsab,expdata%ab,u,intab,a0,b0)


iptr1  = expdata%iptrsalphap(intab,intcd,intef)

call alegendre_tensor_eval(expdata%ncoefsalphap, &
  expdata%coefsalphap,iptr1,a0,b0,c0,d0,e0,f0,u,dmu0,dnu0,apval)


if (ifover .eq. 1) then
apval = apval * dnu
endif


end subroutine

subroutine alegendre_expeval_inverse(expdata,dnu,dmu,t,ainv)
implicit double precision (a-h,o-z)

double precision, intent(in)       :: dnu,dmu,t
type(alegendre_expansion_data)     :: expdata
double precision, intent(out)      :: ainv

ainv  = 0


k           = expdata%k
ifsmalldmu  = expdata%ifsmalldmu
ifover      = expdata%ifover
ifbetas     = expdata%ifbetas
ifinverse   = expdata%ifinverse

if (ifover .eq. 1) then
dnu0   = 1/dnu
else
dnu0   = dnu
endif

call alegendre_evalabc(ifsmalldmu,dnu,dmu,a,b,c)
call compute_dmu0(ifsmalldmu,dnu,dmu,dmu0)
call alegendre_findint(expdata%nintsef,expdata%ef,dnu0,intef,e0,f0)
call alegendre_findint(expdata%nintscd,expdata%cd,dmu0,intcd,c0,d0)

u     = 0
intab = 1
a0    = expdata%ab(1,1)
b0    = expdata%ab(2,1)
iptr  = expdata%iptrsalpha(intab,intcd,intef)
call alegendre_tensor_eval(expdata%ncoefsalpha,expdata%coefsalpha, &
   iptr,a0,b0,c0,d0,e0,f0,u,dmu0,dnu0,aval)

if (ifover .eq. 1) then
aval = aval * dnu 
endif


ainv0  = aval 
binv0  = pi/2*(dnu-dmu) + 2*pi

u = (t-ainv0)/(binv0-ainv0)
call alegendre_findint(expdata%nintsinv,expdata%abinv,u,intinv,a0,b0)
iptr1 = expdata%iptrsalphainv(intinv,intcd,intef)

call alegendre_tensor_eval(expdata%ncoefsalphainv,expdata%coefsalphainv, &
   iptr1,a0,b0,c0,d0,e0,f0,u,dmu0,dnu0,ainv)


if (ifover .eq. 1) then
ainv  = ainv*dnu
endif

end subroutine





subroutine alegendre_expeval00(expdata,dnu,dmu,aval,apval,appval)
implicit double precision (a-h,o-z)

type(alegendre_expansion_data)     :: expdata

!
!  Evaluate alpha and its first two derivatives at the turning point
!

aval   = 0
apval  = 0
bval   = 0


k           = expdata%k
ifsmalldmu  = 0
ifover      = expdata%ifover
ifbetas     = expdata%ifbetas


if (ifover .eq. 1) then
dnu0   = 1/dnu
else
dnu0   = dnu
endif

call alegendre_evalabc(ifsmalldmu,dnu,dmu,a,b,c)
call compute_dmu0(ifsmalldmu,dnu,dmu,dmu0)


u = 0
call alegendre_findint(expdata%nintsef,expdata%ef,dnu0,intef,e0,f0)
call alegendre_findint(expdata%nintscd,expdata%cd,dmu0,intcd,c0,d0)
call alegendre_findint(expdata%nintsab,expdata%ab,u,intab,a0,b0)

iptr1  = expdata%iptrsalpha(intab,intcd,intef)
iptr2  = expdata%iptrsalphap(intab,intcd,intef)

call alegendre_tensor_eval2(expdata%ncoefsalpha,expdata%coefsalpha,expdata%ncoefsalphap, &
  expdata%coefsalphap,iptr1,iptr2,a0,b0,c0,d0,e0,f0,u,dmu0,dnu0,aval,apval)


if (ifover .eq. 1) then
aval  = aval  * dnu
apval = apval * dnu
endif

iptr = expdata%iptrsalphapp(intcd,intef)
call alegendre_tensor_eval2d(expdata%ncoefsalphapp,expdata%coefsalphapp,&
  iptr,c0,d0,e0,f0,dmu0,dnu0,appval)

if (ifover .eq. 1) then
appval = appval * dnu**2
endif


end subroutine


subroutine alegendre_evalabc(ifsmalldmu,dnu,dmu,a,b,c)
implicit double precision (a-h,o-z)


if (ifsmalldmu .eq. 0) then
deta    = sqrt(dmu**2-0.25d0)
dlambda = dnu+0.5d0
a       = asin(deta/dlambda) 
b       = pi/2
c       = asin(deta/dlambda)/100.0d0
else
a       = 1.0d0/dnu**1.1d0
b       = pi/2
endif

end subroutine


subroutine compute_dmu(ifsmalldmu,dnu,dmu,dmu0)
implicit double precision (a-h,o-z)


if (ifsmalldmu .eq. 1 ) then
dmu = dmu0 
else
dmu = 1 + dmu0 * (dnu-1)
endif


end subroutine


subroutine compute_dmu0(ifsmalldmu,dnu,dmu,dmu0)
implicit double precision (a-h,o-z)


if (ifsmalldmu .eq. 1 ) then
dmu0 = dmu
else
dmu0 = (dmu-1)/(dnu-1)
endif

end subroutine



subroutine alegendre_tensor_eval(ncoefs,coefs,iptr, &
  a,b,c,d,e,f,x,y,z,val)
implicit double precision (a-h,o-z)

double precision :: coefs(ncoefs)
double precision :: polsx(0:100),polsy(0:100),polsz(0:100)

nz      = coefs(iptr)
ny      = coefs(iptr+1)
nx      = coefs(iptr+2)

xx      = (x-(b+a)/2)*2/(b-a)
yy      = (y-(d+c)/2)*2/(d-c)
zz      = (z-(f+e)/2)*2/(f-e)


call alegendre_chebs(xx,nx,polsx)
call alegendre_chebs(yy,ny,polsy)
call alegendre_chebs(zz,nz,polsz)

iptr0     = iptr+3
val       = 0

do k=0,nz
n      = coefs(iptr0)
iptr0  = iptr0+1
do j=0,n
m     = coefs(iptr0)
iptr0  = iptr0+1
do i=0,m
val    = val + coefs(iptr0)*polsx(i)*polsy(j)*polsz(k)
iptr0  = iptr0+1
end do
end do
end do

end subroutine



subroutine alegendre_tensor_eval2(ncoefs1,coefs1,ncoefs2,coefs2,iptr1,iptr2, &
  a,b,c,d,e,f,x,y,z,val1,val2)
implicit double precision (a-h,o-z)

double precision :: coefs1(ncoefs1),coefs2(ncoefs2)
double precision :: polsx(0:100),polsy(0:100),polsz(0:100)

nz1     = coefs1(iptr1)
ny1     = coefs1(iptr1+1)
nx1     = coefs1(iptr1+2)

nz2     = coefs2(iptr2)
ny2     = coefs2(iptr2+1)
nx2     = coefs2(iptr2+2)

nx      = max(nx1,nx2)
ny      = max(ny1,ny2)
nz      = max(nz1,nz2)

xx      = (x-(b+a)/2)*2/(b-a)
yy      = (y-(d+c)/2)*2/(d-c)
zz      = (z-(f+e)/2)*2/(f-e)


call alegendre_chebs(xx,nx,polsx)
call alegendre_chebs(yy,ny,polsy)
call alegendre_chebs(zz,nz,polsz)

iptr      = iptr1+3
val1      = 0

do k=0,nz1
n     = coefs1(iptr)
iptr  = iptr+1
do j=0,n
m     = coefs1(iptr)
iptr  = iptr+1
do i=0,m
val1  = val1 + coefs1(iptr)*polsx(i)*polsy(j)*polsz(k)
iptr  = iptr+1
end do
end do
end do

iptr      = iptr2+3
val2      = 0

do k=0,nz2
n     = coefs2(iptr)
iptr  = iptr+1
do j=0,n
m     = coefs2(iptr)
iptr  = iptr+1
do i=0,m
val2  = val2 + coefs2(iptr)*polsx(i)*polsy(j)*polsz(k)
iptr  = iptr+1
end do
end do
end do

end subroutine


subroutine alegendre_tensor_eval2d(ncoefs,coefs,iptr0,a,b,c,d,x,y,val)
implicit double precision (a-h,o-z)

double precision :: coefs(ncoefs)
double precision :: polsx(100),polsy(100)

nx     = coefs(iptr0)
ny     = coefs(iptr0+1)

xx = (x-(b+a)/2)*2/(b-a)
yy = (y-(d+c)/2)*2/(d-c)

call alegendre_chebs(xx,nx,polsx)
call alegendre_chebs(yy,ny,polsy)

val    = 0
iptr   = iptr0+2

do j=0,ny
nj    =  coefs(iptr)
iptr  = iptr+1
do i=0,nj
val     = val  + coefs(iptr)*polsx(i+1)*polsy(j+1)
iptr   = iptr+1
end do
end do


end subroutine




subroutine alegendre_findint(nints,ab,x,int,a,b)
implicit double precision (a-h,o-z)

integer             :: int,nints
double precision    :: ab(2,nints),x,a,b

integer             :: int0

!eps0 = epsilon(0.0d0)

! intl   = 1
! intr   = nints

do int = 1,nints-1
b = ab(2,int)
if (x .le. b) exit
end do

int0 = int
a = ab(1,int)
b = ab(2,int)

end subroutine



subroutine alegendre_chebs(x,n,pols)
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


end module
