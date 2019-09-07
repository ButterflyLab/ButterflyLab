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
!  This file is for the MATLAB interface for James Bremer’s code to evaluate 
!  the Bessel functions of the first and second
!  kinds of real orders and arguments.  For the most part, it uses the precomputed
!  expansions found in the file bessel_data.f90 to do so.
!
!  More specifically, bessel_data.f90 stores representations of three different
!  functions:
!
!  - a nonoscillatory phase function \alpha(\nu, t), defined in the oscillatory region
!      where t > \sqrt(\nu^2-1/4), such that
!
!                             (   2   )^(1/2)   \cos ( \alpha(\nu, t) - pi/2 \nu)
!        J_\nu(t)\sqrt(t)   = ( ----- )         ----------------------------------      (1)
!                             ( \pi t )              \sqrt( \alpha'(\nu,t) ) 
!
!      and
!
!                             (   2   )^(1/2)   \sin ( \alpha(\nu, t) - pi/2 \nu)
!        Y_\nu(t)\sqrt(t)   = ( ----- )         ----------------------------------      (2)
!                             ( \pi t )              \sqrt( \alpha'(\nu,t) ) 
!
!  - the logarithm of J_\nu(t),  defined in the nonoscillatory region 
!    t <= \sqrt(\nu^2-1/4); and
!
!  - the logarithm of -Y_\nu(t), defined in the nonoscillatory region 
!    t <= \sqrt(\nu^2-1/4).
! 
!  These precomputed expansions are supplemented with other approaches in a few cases:
!
!  - In order to evaluate J_\nu(t) and Y_\nu(t) when \nu >= 5 and the argument is large 
!    (t > 1,000,000 \nu ), an asymptotic expansion of the phase function 
!    \alpha(t) is used;
!
!  - In order to evaluate J_\nu(t) and Y_\nu(t) for \nu >= 5 and very small t 
!      (t < \nu / 1,000,000), Debye's asymptotic expansions are used.
!
!  - When \nu <= 5 and t <= 1, series expansions of J_\nu and Y_\nu are used.
!
!
!  The following subroutines should be regarded as user-callable:
!
!    bessel_eval - evaluate the Bessel functions at a specified point (\nu,t)
!      with \nu >= 0 and t > 0.  
!
!      In the event that (\nu,t) is in the oscillatory region (when t >= \sqrt(\nu^2-1/4)), 
!      the values of \alpha(\nu,t) and its derivative with respect to t are also returned.  
!
!      In the case that (\nu,t) is in the nonoscillatory region (when t < \sqrt(\nu^2-1/4)),
!      the values of \log( J_\nu(t)) and \log( - Y_\nu(t)) are also returned.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "fintrf.h"

! standard header for mex functions
subroutine mexfunction(nlhs, plhs, nrhs, prhs)

  ! access bessel module
  use bessel_mod

  implicit double precision (a-h,o-z)

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
  double precision :: dnu
  double precision, allocatable :: ts(:,:),alphas(:,:),alphaders(:,:), & 
       vallogjs(:,:),vallogys(:,:),valjs(:,:),valys(:,:)

  mwPointer    :: t_ptr

  double precision :: dsize


  ! output variable
  double precision :: alpha,alphader,vallogj,vallogy,valj,valy

  character (len=80) :: message

  m     = mxGetM(prhs(2))
  n     = mxGetN(prhs(2))

  ! copy the right-hand side argument in matlab to dnu and t
  allocate(ts(m,n),alphas(m,n),alphaders(m,n),vallogjs(m,n),vallogys(m,n), &
       valjs(m,n),valys(m,n))

  call mxCopyPtrToReal8(mxGetPr(prhs(1)), dnu, 1_mwS)
  call mxCopyPtrToReal8(mxGetPr(prhs(2)), ts, m*n)

  ! initialize the left-hand side of the function
  plhs(1) = mxCreateDoubleMatrix(m, n, 0)
  plhs(2) = mxCreateDoubleMatrix(m, n, 0)
  plhs(3) = mxCreateDoubleMatrix(m, n, 0)
  plhs(4) = mxCreateDoubleMatrix(m, n, 0)
  plhs(5) = mxCreateDoubleMatrix(m, n, 0)
  plhs(6) = mxCreateDoubleMatrix(m, n, 0)

  ! compute the fast Bessel function evaluation 
  do j=1,n
     do i=1,m
        t = ts(i,j)  
        call bessel_eval(dnu,t,alpha,alphader,vallogj,vallogy,valj,valy)
        alphas(i,j)    = alpha
        alphaders(i,j) = alphader
        vallogjs(i,j)  = vallogj
        vallogys(i,j)  = valogy
        valjs(i,j)     = valj
        valys(i,j)     = valy
     end do
  end do


  ! copy the output to the left-hand side
  call mxCopyReal8ToPtr(alphas, mxGetPr(plhs(1)), n*m) 
  call mxCopyReal8ToPtr(alphaders, mxGetPr(plhs(2)), n*m) 
  call mxCopyReal8ToPtr(vallogjs, mxGetPr(plhs(3)), n*m) 
  call mxCopyReal8ToPtr(vallogys,mxGetPr(plhs(4)), n*m) 
  call mxCopyReal8ToPtr(valjs,mxGetPr(plhs(5)), n*m) 
  call mxCopyReal8ToPtr(valys,mxGetPr(plhs(6)), n*m) 


  deallocate(ts,alphas,alphaders,vallogjs,vallogys,valjs,valys)

end subroutine mexfunction
