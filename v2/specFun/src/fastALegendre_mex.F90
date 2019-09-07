!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Copyright 2017 by Haizhao Yang, 2017    
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
!  the *normalized* associated Legendre function of the first kind of 
!  degree nu and order -dmu at the point cos(t).  That is, return the value of
!
!
!     (  (dnu + 1/2) Gamma(dnu+dmu+1) )
!     (  ---------------------------- )^(1/2)  P_dnu^{-dmu}(cos(t)),                      (1)
!     (        Gamma(dnu-dmu+1)       )
!
!  where P_dnu^{-dmu} is the associated Legendre function of the first kind of
!  degree dnu and order -dmu.    The normalization factor in (1) is chosen
!  so as to make the normalized associated Legendre functions of integer orders
!  have unit L^2((-1,1)) norm.
!
!  When (dnu,dmu,t) is in the oscillatory region (that is, when dmu <= 1/2 or
!  when dmu > 1/2 and t < arcsin( sqrt(dmu^2-1/4) / (dnu + 1/2) ) ), this
!  routine also returns the values of a nonoscillatory phase function for
!  the associated Legendre differential equation and its derivative.
!
!  When (dnu,dmu,t) is not in the oscillatory region, it returns the 
!  logarithm of the function (1).
!
!  Input parameters:
!    dnu - the degree of the associated Legendre function to evaluate
!    dmu - the order of the associated Legendre function ot evaluate
!    t - the argument
!
!  Output parameters:
!    alpha - the value of a nonoscillatory phase function
!    alphader - the value of the derivative of the nonoscillatory phase function
!    vallog - the value of the logarithm of (1)
!    val - the value of the normalized associated Legendre function
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "fintrf.h"

! standard header for mex functions
subroutine mexfunction(nlhs, plhs, nrhs, prhs)

  ! access alegendre module
  use alegendre_mod

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
  double precision :: dnu, dmu
  double precision, allocatable :: ts(:,:),alphas(:,:),alphaders(:,:), & 
       vallogps(:,:),valps(:,:),vallogqs(:,:),valqs(:,:)

  mwPointer    :: t_ptr


  ! output variable
  double precision :: alpha,alphader,vallog,val

  character (len=80) :: message

  m     = mxGetM(prhs(3))
  n     = mxGetN(prhs(3))

  ! copy the right-hand side argument in matlab to dnu, dmu and t
  allocate(ts(m,n),alphas(m,n),alphaders(m,n),vallogps(m,n),valps(m,n))
  allocate(vallogqs(m,n),valqs(m,n))

  call mxCopyPtrToReal8(mxGetPr(prhs(1)), dnu, 1_mwS)
  call mxCopyPtrToReal8(mxGetPr(prhs(2)), dmu, 1_mwS)
  call mxCopyPtrToReal8(mxGetPr(prhs(3)), ts, m*n)

  ! initialize the left-hand side of the function
  plhs(1) = mxCreateDoubleMatrix(m, n, 0)
  plhs(2) = mxCreateDoubleMatrix(m, n, 0)
  plhs(3) = mxCreateDoubleMatrix(m, n, 0)
  plhs(4) = mxCreateDoubleMatrix(m, n, 0)
  plhs(5) = mxCreateDoubleMatrix(m, n, 0)
  plhs(6) = mxCreateDoubleMatrix(m, n, 0)

  ! compute the fast *normalized* associated Legendre function 
  do j=1,n
     do i=1,m
        t = ts(i,j)  
        call alegendre_eval(dnu,dmu,t,alpha,alphader,vallogp,vallogq,valp,valq)
        alphas(i,j)    = alpha
        alphaders(i,j) = alphader
        vallogps(i,j)  = vallogp
        valps(i,j)     = valp
        vallogqs(i,j)  = vallogq
        valqs(i,j)     = valq
     end do
  end do


  ! copy the output to the left-hand side
  call mxCopyReal8ToPtr(alphas, mxGetPr(plhs(1)), n*m) 
  call mxCopyReal8ToPtr(alphaders, mxGetPr(plhs(2)), n*m) 
  call mxCopyReal8ToPtr(vallogps, mxGetPr(plhs(3)), n*m)
  call mxCopyReal8ToPtr(vallogqs, mxGetPr(plhs(4)), n*m)
  call mxCopyReal8ToPtr(valps,mxGetPr(plhs(5)), n*m) 
  call mxCopyReal8ToPtr(valqs,mxGetPr(plhs(6)), n*m) 

  deallocate(ts,alphas,alphaders,vallogps,valps,valqs,vallogqs)

end subroutine mexfunction
