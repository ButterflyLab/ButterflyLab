!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  This file contains code for approximating solutions of Jacobi's differential equations
!  using Taylor series expansions. 
!
!  More explicitly, these routines evaluate the functions
!
!   ~(da,db)    (da,db)                    (da,db)
!   P   (t)  = P        ( cos(t) )  r(t)  C                                               (1)
!    dnu        dnu                        dnu
!
!  where  P_dnu^(da,db) is the Jacobi polynomial of the first kind,
!                                                  
!     r(t) = cos(t/2) ^(db+1/2) sin(t/2) ^(da+1/2)
!
!  and  C_dnu^(da,db) is a normalization constant.  The normalization constant is
!  chosen so that the basis of solutions of Jacobi's differential equation consisting
!  of (1) and the corresponding Q function has Wronskian 1.  The L^2 normalized
!  Jacobi polynomials are obtained by multiplying (1) by the constant
!
!       sqrt((2*dnu+da+db+1)/pi) .
!
!  Throughout, we let p denote the quantity  dnu +(da+db+1)/2.  The square of
!  the amplitude function is
!
!
!      (da, db)       ~(da,db)   ^2    ~(da,db)  ^ 2
!     M        (t) =  P      (t)    +  Q     (t)     .                                   (2)
!      dnu             dnu               dnu
!
!  The following subroutines should be regarded as publicly callable:
!
!    jacobi_taylor_p - evaluate (1) and its first two derivatives when dnu is of
!      modest size and t is close to 0
!
!    jacobi_taylor_q - evaluate the q function corresponding to (1) and its first
!      two derivatives when dnu is of modest size and t is close to 0
!
!    jacobi_taylor_m - evaluate the amplitude (2) and its fiest two derivatives
!      when dnu is of modest szie and t is cloe to 0 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine jacobi_taylor_p(dnu,da,db,t,val,der,der2)
implicit double precision (a-h,o-z)

!  
!  Evaluate the function (1) at a point near 0 when dnu is of modest size.
!

data pi / 3.14159265358979323846264338327950288d0 /

rval  = cos(t/2)**(db+0.5d0) * sin(t/2) ** (da+0.5d0)
rder  = (Cos(t/2.)**(-0.5d0 + db)*(da - db + (1 + da + db)*Cos(t))*Sin(t/2.)**(-0.5d0 + da))/4d0
rder2 =         (Cos(t/2)**(-1.5d0 + db)*(-3 - 2*da + 3*da**2 - 2*db - 2*da*db + 3*db**2 + &
           4*(da**2 - db**2)*Cos(t) + (1 + da + db)**2*Cos(2*t))*Sin(t/2)**(-1.5d0 + da))/32.0d0

x = cos(t)

call jacobi_taylor_j1(dnu,da,db,x,fval,fder,fder2)

gval  = fval
gder  = -sin(t)*fder
gder2 = -cos(t)*fder + sin(t)**2*fder2

fval  = gval*rval
fder  = gder*rval + gval*rder
fder2 = gder2*rval + 2*gder*rder + gval*rder2

!
!  Now scale by the Wronskian
!

call jacobi_gammaratio1(dnu+1,da,dd1)
call jacobi_gammaratio1(dnu+db+1,da,dd2)
dconst = dd1/dd2 * 1/pi
dconst = 1/sqrt(dconst)

val  = fval*dconst
der  = fder*dconst
der2 = fder2*dconst

end subroutine



subroutine jacobi_taylor_q(dnu,da,db,t,val,der,der2)
implicit double precision (a-h,o-z)

!  
!  Evaluate the scaled, Wronskain normalized Jacobi function of the
!  second kind at a point t near 0.
!

data pi / 3.14159265358979323846264338327950288d0 /

rval  = cos(t/2)**(db+0.5d0) * sin(t/2) ** (da+0.5d0)
rder  = (Cos(t/2.)**(-0.5d0 + db)*(da - db + (1 + da + db)*Cos(t))*Sin(t/2.)**(-0.5d0 + da))/4d0
rder2 =         (Cos(t/2)**(-1.5d0 + db)*(-3 - 2*da + 3*da**2 - 2*db - 2*da*db + 3*db**2 + &
           4*(da**2 - db**2)*Cos(t) + (1 + da + db)**2*Cos(2*t))*Sin(t/2)**(-1.5d0 + da))/32.0d0

x = cos(t)

call jacobi_taylor_j4(dnu,da,db,x,fval,fder,fder2)

gval  = fval
gder  = -sin(t)*fder
gder2 = -cos(t)*fder + sin(t)**2*fder2

fval  = gval*rval
fder  = gder*rval + gval*rder
fder2 = gder2*rval + 2*gder*rder + gval*rder2

!
!  Now scale by the Wronskian
!

call jacobi_gammaratio1(dnu+1,da,dd1)
call jacobi_gammaratio1(dnu+db+1,da,dd2)
dconst = dd1/dd2 * 1/pi
dconst = 1/sqrt(dconst)

val  = fval*dconst
der  = fder*dconst
der2 = fder2*dconst

end subroutine



subroutine jacobi_taylor_j1(dnu,da,db,x,val,der,der2)
implicit double precision (a-h,o-z)

!
!  Evaluate the Jacobi function J1 at a point  x  close to 1.
!

data pi / 3.14159265358979323846264338327950288d0 /

!
!  Evaluate the hypergeometric series
!

maxterms = 100
eps0     = 1.0d-30
dterm    = 1.0d0
xx       = (x-1)/(x+1)

val0      = 0 
der0      = 0
der20     = 0

do k=1,maxterms

val0   = val0  + dterm
der0   = der0  + dterm*2*(k-1)/(x**2-1)
der20  = der20 + dterm*4*(k-1)*(k-1-x)/(x**2-1)**2

dterm = dterm * xx/k
dterm = dterm * (k-dnu-db-1)*(k-dnu-1)/(k+da) 

!if (abs(dterm) .lt. eps0*abs(val0)) exit

end do

!
!  Mulitply by ((1+x)/2)**n
!


val  = val0  *((1+x)/2)**dnu
der  = der0 * ((1+x)/2)**dnu + val0 * 2.0d0**(-dnu)*dnu*(1+x)**(dnu-1)
der2 = ((1 + x)**(-2 + dnu)*((-1 + dnu)*dnu*val0 + &
        (1 + x)*(2*der0*dnu + der20*(1 + x))))/2**dnu

!
!  Scale by the appropriate  constant
!

call jacobi_gammaratio1(dnu+1,da,dd1)
dconst = dd1 / gamma(1+da)

val  = val * dconst
der  = der * dconst
der2 = der2 * dconst

end subroutine


subroutine jacobi_taylor_j2(dnu,da,db,x,val,der,der2)
implicit double precision (a-h,o-z)

!  
!  Evaluate the Jacobi function J2 at a point close to 1.
!

data pi / 3.14159265358979323846264338327950288d0 /

maxterms = 100
eps0     = 1.0d-30

!
!  Evaluate the hypergeometric series 2F1(-n-a,-n-a-b,1-a,(x-1)/(x+1))
!

dterm    = 1.0d0
xx       = (x-1)/(x+1)
val0     = 0 
der0     = 0
der20    = 0
do k=1,maxterms
val0   = val0  + dterm
der0   = der0  + dterm*2*(k-1)/(x**2-1)
der20  = der20 + dterm*4*(k-1)*(k-1-x)/(x**2-1)**2

dterm = dterm * xx/k
dterm = dterm * (-1-da+k-dnu)*(-1-da-db+k-dnu)/(k-da)
if (abs(dterm) .lt. eps0*abs(val0)) exit
end do

!
!  Multiply by ((1+x)/2)^n
!

val2  = val0  *((1+x)/2)**dnu
der2  = der0 * ((1+x)/2)**dnu + val0 * 2.0d0**(-dnu)*dnu*(1+x)**(dnu-1)
der22 = ((1 + x)**(-2 + dnu)*((-1 + dnu)*dnu*val0 + &
        (1 + x)*(2*der0*dnu + der20*(1 + x))))/2**dnu

!
!  Multiply by ((1+x)/(1-x))^da
!

val3  = val2 * ((1+x)/(1-x))**da
der3  = (((1 + x)/(1 - x))**da*(-2*da*val2 + der2*(-1 + x**2)))/(-1 + x**2)
der32 =         (((1 + x)/(1 - x))**da*(4*da*val2*(da + x) +  &
           (-1 + x**2)*(-4*da*der2 + der22*(-1 + x**2))))/(-1 + x**2)**2

!
!  Scale by the appropriate  constant
!

call jacobi_gammaratio1(dnu+db+1,da,dd1)
dconst = 1/dd1 *1/ gamma(1-da)


val  = val3  * dconst
der  = der3  * dconst
der2 = der32 * dconst

end subroutine



recursive subroutine jacobi_taylor_j4(dnu,da,db,x,val,der,der2)
implicit double precision (a-h,o-z)

!  
!  Evaluate the appropriate function of the second kind at the point x.
!

double precision :: whts(16),vals(16),ders(16),der2s(16),xs(16)

data pi / 3.14159265358979323846264338327950288d0 /

data ncheb / 16 /

data xs / -0.100000000000000000000000000000000000D+01, &
          -0.978147600733805637928566747869599532D+00, &
          -0.913545457642600895502127571985317191D+00, &
          -0.809016994374947424102293417182819080D+00, &
          -0.669130606358858213826273330686780548D+00, &
          -0.500000000000000000000000000000000096D+00, &
          -0.309016994374947424102293417182818984D+00, &
          -0.104528463267653471399834154802498051D+00, &
           0.104528463267653471399834154802498147D+00, &
           0.309016994374947424102293417182819032D+00, &
           0.499999999999999999999999999999999952D+00, &
           0.669130606358858213826273330686780452D+00, &
           0.809016994374947424102293417182819080D+00, &
           0.913545457642600895502127571985317191D+00, &
           0.978147600733805637928566747869599532D+00, &
           0.100000000000000000000000000000000000D+01  /

! data xs / -0.1000000000000000D+01, &
!           -0.9781476007338056D+00, &
!           -0.9135454576426008D+00, &
!           -0.8090169943749473D+00, &
!           -0.6691306063588579D+00, &
!           -0.4999999999999998D+00, &
!           -0.3090169943749473D+00, &
!           -0.1045284632676533D+00, &
!            0.1045284632676537D+00, &
!            0.3090169943749475D+00, &
!            0.5000000000000001D+00, &
!            0.6691306063588582D+00, &
!            0.8090169943749475D+00, &
!            0.9135454576426009D+00, &
!            0.9781476007338057D+00, &
!            0.1000000000000000D+01  /

!
!  When da of small magnitude, use interpolation to evaluate j4
!

if (abs(da) .lt. 0.01d0) then
eps0  = epsilon(0.0d0)
a     = -0.097d0
b     =  0.097d0
xx    = (2*da - (b+a) ) /(b-a)

do i=1,ncheb
alpha = xs(i)*(b-a)/2+ (b+a)/2
call jacobi_taylor_j4(dnu,alpha,db,x,vals(i),ders(i),der2s(i))
end do

!
!  do the interpolation
!

dsign = 1.0d0
do i=1,ncheb
diff    = xx - xs(i)
if (abs(diff) .lt. eps0) then
val  = vals(i)
der  = ders(i)
der2 = der2s(i)
return
endif

whts(i) = dsign/diff
dsign   = -dsign
end do
whts(1)     = whts(1)/2 
whts(ncheb) = whts(ncheb)/2

dd = sum(whts)
val  = sum(vals*whts)/dd
der  = sum(ders*whts)/dd
der2 = sum(der2s*whts)/dd

else

dd1 = gamma(1-da)/(gamma(0.5d0+da)*gamma(0.5d0-da))
dd2 = gamma(1-da)/pi

call jacobi_taylor_j1(dnu,da,db,x,val1,der11,der12)
call jacobi_taylor_j2(dnu,da,db,x,val2,der21,der22)

val  = dd1*val1-dd2*val2
der  = dd1*der11-dd2*der21
der2 = dd1*der12-dd2*der22

val  = val * gamma(da)
der  = der * gamma(da)
der2 = der2* gamma(da)

endif

end subroutine






subroutine jacobi_taylor_m(dnu,da,db,t,val,der,der2)
implicit double precision (a-h,o-z)

!  
!  Evaluate the scaled, Wronskain normalized Jacobi function of the
!  second kind at a point t near 0.
!

data pi / 3.14159265358979323846264338327950288d0 /

call jacobi_taylor_q(dnu,da,db,t,qval,qder,qder2)
call jacobi_taylor_p(dnu,da,db,t,pval,pder,pder2)

val  = pval**2+qval**2
der  = 2*pval*pder + 2*qval*qder
der2 = 2*pder*pder + 2*pval*pder2 + 2*qder*qder+2*qval*qder2

end subroutine

