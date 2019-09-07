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
!
    !
    !  Evaluate only the derivative of alpha, and only for dnu > 10 and in
    !  the oscillatory regime.
    !  
    !
!
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
  mwSize       :: ms 
  integer,  parameter :: mwS = kind(ms)

  ! input variable
  double precision :: dnu, dmu

  ! output variable
  double precision :: alphader

  character*120 line
  double precision :: t

  ! copy the right-hand side argument in matlab to dnu, dmu
  call mxCopyPtrToReal8(mxGetPr(prhs(1)), dnu, 1_mwS)
  call mxCopyPtrToReal8(mxGetPr(prhs(2)), dmu, 1_mwS)
  call mxCopyPtrToReal8(mxGetPr(prhs(3)), t, 1_mwS)

  ! compute the smooth phase and its first two derivatives
  call  alegendre_eval00(dnu,dmu,t,alphader)

  ! initialize the left-hand side of the function as a 1-by-1 matrix
  plhs(1) = mxCreateDoubleMatrix(1_mwS, 1_mwS, 0)

  ! copy the output to the left-hand side  
  call mxCopyReal8ToPtr(alphader,mxGetPr(plhs(1)),1_mwS) 

end subroutine mexfunction
