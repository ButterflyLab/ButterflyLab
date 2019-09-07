!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  This module contains miscellaneous utility routines which can be divided into the 
!  following categories:
!
!  Printing routines:
!  -----------------
!
!  The folllowing routines output the contents of arrays.  They have two principal
!  advantages over simply using "print": (1) the output of these routines is also
!  written into the file fort.13 for later consideration, (2) the output is formatted 
!  in a consistent fashion.
!
!    prin2 - output an array of doubles with 7 digits of each element displayed 
!    prind - output an array of doubles with 15 digits of each element displayed
!    prini - output an array of integers
!    prinl - output an array of long integers
!    prina - output a string
!
!  Plotting:
!  ------------------
!
!  This module contains the following very simple plotting routines:
!
!    gnuplot_points - output a GNUplot script and data file which generates
!     a simple two-dimensional plot of a collection of points
!
!  Other routines:
!  ---------------
!
!    elapsed - return the wall clock time in seconds which has elapsed since some
!     arbitrary point in the past 
!
!    mach_zero - return machine zero
!
!    insort0 - sort an array of real numbers
!
!    qrsolv - solve a system of linear equations via a QR decomposition
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


module utils

interface prin2
   module procedure prin2_0
   module procedure prin2_1
   module procedure prin2_2
end interface prin2

interface prind
   module procedure prind_0
   module procedure prind_1
   module procedure prind_2
end interface prind

interface prini
   module procedure prini_0
   module procedure prini_1
   module procedure prini_2
end interface prini

interface prinl
   module procedure prinl_0
   module procedure prinl_1
   module procedure prinl_2
end interface prinl

contains

subroutine prin2_0(str,a)
implicit double precision (a-h,o-z)

double precision a
character(len=*), intent(in) :: str

print *,str
print "(8(2x,e15.7))",a

write (13,*) str
write (13,"(8(2x,e15.7))") a

end subroutine

subroutine prin2_1(str,a)
implicit double precision (a-h,o-z)

double precision, intent (in) :: a(:)
character(len=*), intent(in) :: str

print *,str
print "(8(2x,e15.7))",a

write (13,*) str
write (13,"(8(2x,e15.7))") a

end subroutine

subroutine prin2_2(str,a)
implicit double precision (a-h,o-z)

double precision, intent (in) :: a(:,:)
character(len=*), intent(in) :: str

print *,str
print "(8(2x,e15.7))",a

write (13,*) str
write (13,"(8(2x,e15.7))") a

end subroutine

subroutine prind_0(str,a)
implicit double precision (a-h,o-z)

double precision :: a
character(len=*), intent(in) :: str

print *,str
print "(5(2x,e24.16))",a

write (13,*) str
write (13,"(5(2x,e24.16))") a

end subroutine

subroutine prind_1(str,a)
implicit double precision (a-h,o-z)

double precision,intent(in) :: a(:)
character(len=*), intent(in) :: str

print *,str
print "(5(2x,e24.16))",a

write (13,*) str
write (13,"(5(2x,e24.16))") a

end subroutine


subroutine prind_2(str,a)
implicit double precision (a-h,o-z)

double precision,intent(in) :: a(:,:)
character(len=*), intent(in) :: str

print *,str
print "(5(2x,e24.16))",a

write (13,*) str
write (13,"(5(2x,e24.16))") a

end subroutine



subroutine prini_0(str,a)
implicit double precision (a-h,o-z)

integer,intent(in) :: a
character(len=*), intent(in) :: str

print *,str
print "(8(2x,I8))",a

write (13,*) str
write (13,"(8(2x,I8))") a

end subroutine

subroutine prini_1(str,a)
implicit double precision (a-h,o-z)

integer,intent(in) :: a(:)
character(len=*), intent(in) :: str

print *,str
print "(8(2x,I8))",a

write (13,*) str
write (13,"(8(2x,I8))") a

end subroutine



subroutine prini_2(str,a)
implicit double precision (a-h,o-z)

integer,intent(in) :: a(:,:)
character(len=*), intent(in) :: str

print *,str
print "(8(2x,I8))",a

write (13,*) str
write (13,"(8(2x,I8))") a

end subroutine

subroutine prina(str)
implicit double precision (a-h,o-z)

character(len=*), intent(in) :: str

print *,str
write (13,*) str

end subroutine


subroutine prinl_0(str,a)
implicit double precision (a-h,o-z)

integer,intent(in) :: a
character(len=*), intent(in) :: str

print *,str
print "(6(2x,I13))",a

write (13,*) str
write (13,"(6(2x,I13))") a

end subroutine

subroutine prinl_1(str,a)
implicit double precision (a-h,o-z)

integer,intent(in) :: a(:)
character(len=*), intent(in) :: str

print *,str
print "(6(2x,I13))",a

write (13,*) str
write (13,"(6(2x,I13))") a

end subroutine

subroutine prinl_2(str,a)
implicit double precision (a-h,o-z)

integer,intent(in) :: a(:,:)
character(len=*), intent(in) :: str

print *,str
print "(6(2x,I13))",a

write (13,*) str
write (13,"(6(2x,I13))") a

end subroutine




subroutine mach_zero(zero_mach)
implicit double precision (a-h,o-z)
data eps0 / -1.0d0 /
save

!
!  Only perform the computation once and save the results.
!
if (eps0 .ne. -1.0d0) then
   zero_mach=eps0
   return
endif
!
!  Approximate machine zero in some reasonable way.
!
zero_mach=1
d = 1

do while ( (d + zero_mach/2) .ne. d )
   zero_mach=zero_mach/2
end do

end subroutine



subroutine elapsed(t)
implicit double precision (a-h,o-z)
integer*8 i,irate
real t1
call system_clock(i,irate)

dd = i
dd = dd/irate
t = dd
return
end subroutine


subroutine gnuplot_points(title,iplot,n,xs,ys)
implicit double precision (a-h,o-z)
character*1 title(1),AST
dimension xs(n),ys(n)
double precision, allocatable :: zs(:,:)

character*5 filename1
character*9 filename2
character*80 command
character*1 quote

data quote /'"'/, AST / '*' /


 1000 format("gn",I3.3)
 1100 format("gn",I3.3,".dat")
 1200 format((E24.16,2X,E24.16))
 1300 format(" plot [",E24.16,":",E24.16,"] [",E24.16,":",E24.16," ] ",   &
             '"',A9,'"'," with points pointtype 7 pointsize .25")
 1400 format("gnuplot -persist ",A9)
 1500 format(" set title ",'"',80A1,'"')

!
!       Find the length of the title.
!
len = 0
do j=1,10000
if (title(j) .eq. AST) exit
len=len+1
end do

! 
!        Open the data file, gn???.dat, and the GNUplot script file,
!        gn???.
! 
write(filename1,1000) iplot
write(filename2,1100) iplot

iw1=20
iw2=21

open(iw1,FILE=filename1)
open(iw2,FILE=filename2)

!
!       Write the coordinates of the points to the data file and 
!       close it.

write(iw2,1200) (xs(j),ys(j),j=1,n)
close(iw2)

! 
!        Decide on the length of the axes.
! 
xmin = 1d50
xmax =-1d50
ymin = 1d50
ymax =-1d50

! 
do j=1,n
xmax = max(xs(j),xmax)
xmin = min(xs(j),xmin)        
ymax = max(ys(j),ymax)       
ymin = min(ys(j),ymin)
end do
! 
write(iw1,*)    "set terminal qt"
write(iw1,*)    "unset key"
! 
if (len .gt. 0) then
write (iw1,1500) (title(j),j=1,len)
endif

write(iw1,1300) xmin,xmax,ymin,ymax,filename2

close(iw1)

end subroutine


subroutine insort0(k,a)
implicit double precision (a-h,o-z)
integer, intent(in)              :: k
double precision, intent (inout) :: a(k)

if (k .le. 1) return

do i=2,k
val=a(i)
j=i-1
do while (j .ge. 1 .AND. a(j) .gt. val) 
a(j+1)=a(j)
j=j-1
end do
a(j+1)=val
end do
end subroutine


subroutine insort00(k,ia)
implicit double precision (a-h,o-z)
integer, intent(in)              :: k
integer, intent (inout)          :: ia(k)

if (k .le. 1) return

do i=2,k
ival=ia(i)
j=i-1
do while (j .ge. 1 .AND. ia(j) .gt. ival) 
ia(j+1)=ia(j)
j=j-1
end do
ia(j+1)=ival
end do
end subroutine


subroutine qrsolv(a,n,rhs)
implicit double precision (a-h,o-z)

integer            :: n
double precision   :: a(n,n),rhs(n)

!
!  This subroutine uses a version of QR-decomposition to solve the equation
!  A x = b.  Both the input matrix a and the right-hand side are destroyed
!  b ythe routine.
!
!  Input parameters:
!    a - the (n,n) matrix of coefficients
!    n - an integer specifying the size of the system of 
!    rhs - a vector of length n speciying the rhs of the system
!
!  Output parameters:
!
!   rhs - upon return, the solution of the linear system


double precision :: aa(2),u(2,2)

! 
! transpose the input matrix a 
!

size22=0
do i=1,n
do j=1,i
d=a(j,i)
a(j,i)=a(i,j)
a(i,j)=d
size22=size22+a(j,i)**2
size22=size22+a(i,j)**2
end do
end do

! 
! eliminate the upper right triangle
! 
do i=1,n-1
do j=n,i+1,-1
aa(1)=a(i,j-1)
aa(2)=a(i,j)

! call qrrotfnd2(aa,u,size22)

u22=-aa(1)
u12=aa(2)
d=u22**2+u12**2
if(d .lt. size22*1.0d-66) then
u(2,2)=1
u(1,2)=0
u(1,1)=1
u(2,1)=0
else
d=sqrt(d)
u(2,2)=u22/d
u(1,2)=u12/d
u(1,1)=-u(2,2)
u(2,1)=u(1,2)
endif

!call qrrotate(u,a(1,j-1),a(1,j),n,i)

do ii=i,n
d1=u(1,1)*a(ii,j-1)+u(1,2)*a(ii,j)
d2=u(2,1)*a(ii,j-1)+u(2,2)*a(ii,j) 
a(ii,j-1)=d1
a(ii,j)=d2
end do

d1=u(1,1)*rhs(j-1)+u(1,2)*rhs(j)
d2=u(2,1)*rhs(j-1)+u(2,2)*rhs(j)
rhs(j-1)=d1
rhs(j)=d2
end do
end do

!
!  Apply the inverse of the triangular matrix
! 

!call qrtrin(a,rhs,n)

rhs(n)=rhs(n)/a(n,n)
do i=n-1,1,-1
d=0
do j=n,i+1,-1
d=d+a(j,i)*rhs(j)
end do
rhs(i)=(rhs(i)-d)/a(i,i)
end do

return
end subroutine

end module
