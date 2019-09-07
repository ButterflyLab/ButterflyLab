
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  A MATLAB interface for James Bremer's Fortran code.
!
!  By Haizhao Yang, 2017    
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "fintrf.h"

! standard header for mex functions
subroutine mexfunction(nlhs, plhs, nrhs, prhs)

  ! access alegendre module
  use alegendre_mod

  ! number of input arguments, number of output arguments
  integer      :: nlhs, nrhs  
  ! pointer to inputs and outputs
  mwPointer    :: plhs(*), prhs(*) 
  ! get some of the matlab mex functions
  mwPointer    :: mxGetPr, mxCreateDoubleMatrix 
  ! define a size integer so that we can get its type
  mwSize       :: ms, m, n 
  integer,  parameter :: mwS = kind(ms)

  ! input variable
  double precision :: dmu
  integer :: dnu, j, icnt, jcnt

  ! output variable
  double precision, allocatable :: t(:,:), jth(:,:), wht(:,:)
  character*120 line
  double precision :: dsize

  m     = mxGetM(prhs(3))
  n     = mxGetN(prhs(3))

  allocate(t(m,n),jth(m,n),wht(m,n))

  ! copy the right-hand side argument in matlab to dnu, dmu and t
  call mxCopyPtrToReal8(mxGetPr(prhs(1)), dmu, 1_mwS)
  call mxCopyPtrToReal8(mxGetPr(prhs(2)), dnu, 1_mwS)
  call mxCopyPtrToReal8(mxGetPr(prhs(3)), jth, m*n)

  ! compute the number of roots of the associated Legendre function  
  do jcnt=1,n
     do icnt=1,m
        j=INT(jth(icnt,jcnt))
        call  alegendre_jacobi(dmu,dnu,j,t(icnt,jcnt),wht(icnt,jcnt))
     end do
  end do

  ! initialize the left-hand side of the function as a 1-by-1 matrix
  plhs(1) = mxCreateDoubleMatrix(m, n, 0)
  plhs(2) = mxCreateDoubleMatrix(m, n, 0)

  ! copy the output to the left-hand side
  call mxCopyReal8ToPtr(t, mxGetPr(plhs(1)), m*n) 
  call mxCopyReal8ToPtr(wht, mxGetPr(plhs(2)), m*n) 

end subroutine mexfunction
