
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
!  Return the j^{th} largest root of the function 
!
!     (  (dnu + 1/2) Gamma(dnu+dmu+1) )
!     (  ---------------------------- )^(1/2)  P_dnu^{-dmu}(\cos(t)) \sqrt{\sin(t)}       (1)
!     (        Gamma(dnu-dmu+1)       )
!  on the interval (0,\pi).
!
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
  double precision :: dnu,dmu
  integer :: j, icnt, jcnt

  ! output variable
  double precision, allocatable :: t(:,:), jth(:,:)
  character*120 line
  double precision :: dsize

  m     = mxGetM(prhs(3))
  n     = mxGetN(prhs(3))

  allocate(t(m,n),jth(m,n))

  ! copy the right-hand side argument in matlab to dnu, dmu and t
  call mxCopyPtrToReal8(mxGetPr(prhs(1)), dnu, 1_mwS)
  call mxCopyPtrToReal8(mxGetPr(prhs(2)), dmu, 1_mwS)
  call mxCopyPtrToReal8(mxGetPr(prhs(3)), jth, m*n)

  ! compute the number of roots of the associated Legendre function  
  do jcnt=1,n
     do icnt=1,m
        j=INT(jth(icnt,jcnt))
        call  alegendre_proot(dnu,dmu,j,t(icnt,jcnt))
     end do
  end do


  ! write(line,*) 'proot = ', t
  ! k=mexPrintf(line//achar(13))

  ! initialize the left-hand side of the function as a 1-by-1 matrix
  plhs(1) = mxCreateDoubleMatrix(m, n, 0)

  ! copy the output to the left-hand side
  call mxCopyReal8ToPtr(t, mxGetPr(plhs(1)), m*n) 

end subroutine mexfunction
