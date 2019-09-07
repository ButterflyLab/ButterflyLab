!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  A MATLAB interface for James Bremer's Fortran code.
!
!  By Haizhao Yang, 2018    
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
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "fintrf.h"

! standard header for mex functions
subroutine mexfunction(nlhs, plhs, nrhs, prhs)

  ! access alegendre module
  use legendre

  ! number of input arguments, number of output arguments
  integer      :: nlhs, nrhs  
  ! pointer to inputs and outputs
  mwPointer    :: plhs(*), prhs(*) 
  ! get some of the matlab mex functions
  mwPointer    :: mxGetPr, mxCreateDoubleMatrix 
  ! define a size integer so that we can get its type
  mwSize       :: ms 
  integer,  parameter :: mwS = kind(ms)

  ! input variable
  integer :: n

  ! output variable
  double precision, allocatable :: roots(:), weights(:)

  character*120 line
  double precision :: dsize

  ! copy the right-hand side argument in matlab to dnu, dmu
  call mxCopyPtrToReal8(mxGetPr(prhs(1)), dsize, 1)
  n=INT(dsize)

  ! compute the smooth phase and its first two derivatives
  call  legequad(n,roots,weights)

  ! initialize the left-hand side of the function as a 1-by-1 matrix
  plhs(1) = mxCreateDoubleMatrix(1, n, 1)
  plhs(2) = mxCreateDoubleMatrix(1, n, 1)

  ! copy the output to the left-hand side
  call mxCopyReal8ToPtr(roots,mxGetPr(plhs(1)),n)  
  call mxCopyReal8ToPtr(weights,mxGetPr(plhs(2)),n)  

end subroutine mexfunction
