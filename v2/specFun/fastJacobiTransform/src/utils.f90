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
!    prind - output an array of doubles with 15 digits of each element displayed1
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
!      a simple two-dimensional plot of a collection of points
!    psplot_rectangles - output a postscript file which draws a collection
!      of rectangles
!    mathematica_array - generate a text file containing mathematica code
!      which constructs an array
!
!  Other routines:
!  ---------------
!
!    elapsed - return the wall clock time in seconds which has elapsed since some
!     arbitrary point in the past 
!    insort - sort an array of real numbers
!    insorti  - sort an array of integers
!    quicksort - sort an array of real numbers
!    quicksorti - sort an array of integers
!    qrsolv - solve a system of linear equations via a QR decomposition
!    eye - return a (k,k) identity matrix
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


module utils

interface prin2
   module procedure prin2_0
   module procedure prin2_1
   module procedure prin2_2
   module procedure prin2_3
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


subroutine prin2_3(str,a)
implicit double precision (a-h,o-z)

double precision, intent (in) :: a(:,:,:)
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


subroutine psplot_rectangles(filename,nrects,rects)
implicit double precision (a-h,o-z)

character(len=*)               :: filename
integer                        :: nrects
double precision               :: rects(4,nrects)


xmin = 1d50
xmax =-1d50
ymin = 1d50
ymax =-1d50

do i=1,nrects
x1 = rects(1,i)
x2 = rects(2,i)
y1 = rects(3,i)
y2 = rects(4,i)
xmin = min(xmin,x1)
xmin = min(xmin,x2)
xmax = max(xmax,x1)
xmax = max(xmax,x2)
ymin = min(ymin,y1)
ymin = min(ymin,y2)
ymax = max(ymax,y1)
ymax = max(ymax,y2)
end do

xscale = (xmax-xmin)


iw = 66678
open(iw,FILE=filename)

do i=1,nrects

x1 = rects(1,i)
x2 = rects(2,i)
y1 = rects(3,i)
y2 = rects(4,i)

x1 = (x1 - xmin) * 612/(xmax-xmin)
x2 = (x2 - xmin) * 612/(xmax-xmin)
y1 = (y1 - ymin) * 612/(ymax-ymin)
y2 = (y2 - ymin) * 612/(ymax-ymin)

write (iw,*) "newpath"
write (iw,*) x1," ",y1," moveto"
write (iw,*) x2," ",y1," lineto"
write (iw,*) x2," ",y2," lineto"
write (iw,*) x1," ",y2," lineto"
write (iw,*) x1," ",y1," lineto"
write (iw,*) "stroke"

end do


close(iw)

end subroutine


subroutine gnuplot_points(title,iplot,n,xs,ys)
implicit double precision (a-h,o-z)
character*1 title(1),AST
dimension xs(n),ys(n)
double precision, allocatable :: zs(:,:)

character*6 filename1
character*10 filename2
character*80 command
character*1 quote

data quote /'"'/, AST / '*' /


 1000 format("gn",I4.4)
 1100 format("gn",I4.4,".dat")
 1200 format((E24.16,2X,E24.16))
 1300 format(" plot [",E24.16,":",E24.16,"] [",E24.16,":",E24.16," ] ",   &
             '"',A10,'"'," with points pointtype 7 pointsize .25")
 1400 format("gnuplot -persist ",A10)
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


subroutine insort(k,a)
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


subroutine insorti(k,ia)
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


subroutine insorti2(k,ia,idxs)
implicit double precision (a-h,o-z)
integer, intent(in)              :: k
integer, intent (inout)          :: ia(k),idxs(k)

if (k .le. 1) return

do i=2,k
ival   = ia(i)
idxval = idxs(i)
j=i-1
do while (j .ge. 1 .AND. ia(j) .gt. ival) 
ia(j+1)   = ia(j)
idxs(j+1) = idxs(j)
j=j-1
end do
ia(j+1)   = ival
idxs(j+1) = idxval
end do
end subroutine


subroutine quicksort(n,vals)
implicit double precision (a-h,o-z)
dimension istack(2,20000)
dimension vals(1),idxs(1)
!
!       Sort a list of double precision numbers.
!
        if (n .lt. 100) then
        call insort(n,vals)
        return
        endif
!
        maxstack = 10000
        k        = 60
!
        m = 1
        istack(1,1) = 1
        istack(2,1) = n
!
 1000 continue
        if (m .eq. 0) goto 1100
        i1 = istack(1,m)
        i2 = istack(2,m)
        m=m-1
!
        l = i2-i1+1
        if (l .le. k) then
        call insort(l,vals(i1))
        goto 1000
        endif
!
!       Otherwise perform quicksort.
!
        call quicksort01(vals,i1,i2,i3)
!
!       This should never happen, but just in case ...
!
        if (m+2 .ge. maxstack) then
        print *,"quicksort out of memory"
        stop
        endif
!
!       Make sure the smaller half is processed first to reduce storage
!       to O(logn).
!             
        n1 = i3-i1+1
        n2 = i2-i3
!
        if (n2 .lt. n1) then
!
        m = m+1
        istack(1,m) = i1
        istack(2,m) = i3
!
        m = m+1
        istack(1,m) = i3+1
        istack(2,m) = i2
!
        else
!
        m = m+1
        istack(1,m) = i3+1
        istack(2,m) = i2
!
        m = m+1
        istack(1,m) = i1
        istack(2,m) = i3
!
        endif
!
        goto 1000
 1100 continue 
end subroutine


        subroutine quicksort01(vals,i1,i2,i3)
        implicit double precision (a-h,o-z)
        dimension vals(1)
!
!       Randomly choose a pivot index.
!
!        call corrand3(1,r)
        call random_number(r)
        ipiv = i1+floor((i2-i1)*r)
!
!        ipiv = i1+(i2-i1)/2
!
        val  = vals(ipiv)
!
!       Swap the pivot element and the last element.
!
        vals(ipiv) = vals(i2)
        vals(i2)   = val
!
       i3 = i1
!
        do 1000 i=i1,i2-1
!
        if( vals(i) .lt. val) then
        d  = vals(i)
!
        vals(i)  = vals(i3)
        vals(i3) = d       
!
        i3=i3+1
        endif
!
 1000 continue
!
        dd = vals(i3)
        vals(i3) = vals(i2)
        vals(i2) = dd
!
        end subroutine


        subroutine quicksorti(n,ivals)
        implicit double precision (a-h,o-z)
        dimension istack(2,20000),ivals(1)
!
!       Sort a list of integers.
!
        if (n .lt. 60) then
        call insorti(n,ivals)
        return
        endif
!
        maxstack = 10000
        k        = 60
!
        nstack      = 1
        istack(1,1) = 1
        istack(2,1) = n
!
 1000 continue
        if (nstack .eq. 0) goto 1100
        i1 = istack(1,nstack)
        i2 = istack(2,nstack)
        nstack=nstack-1
!
!
        l = i2-i1+1
        if (l .le. k) then
        call insorti(l,ivals(i1))
        goto 1000
        endif
!
!       Otherwise perform quicksort.
!
        call quicksorti0(ivals,i1,i2,i3)
!
!       This should never happen, but just in case ...
!
        if (nstack+2 .ge. maxstack) then
        print *,"quicksort out of memory"
        stop
        endif
!
!       Make sure the smaller half is processed first to reduce storage
!       to O(logn).
!             
        n1 = i3-i1+1
        n2 = i2-(i3+1)+1
!
        if (n2 .lt. n1) then
!
        nstack = nstack+1
        istack(1,nstack) = i1
        istack(2,nstack) = i3
!
        nstack = nstack+1
        istack(1,nstack) = i3+1
        istack(2,nstack) = i2
!
        else
!
        nstack=nstack+1
        istack(1,nstack) = i3+1
        istack(2,nstack) = i2
!
        nstack=nstack+1
        istack(1,nstack) = i1
        istack(2,nstack) = i3
!
        endif
!
        goto 1000
 1100 continue
        end subroutine


        subroutine quicksorti0(ivals,i1,i2,i3)
        implicit double precision (a-h,o-z)
        dimension ivals(1)
!
!       Randomly choose a pivot index.
!
!        call corrand3(1,r)
        call random_number(r)
        ipiv = i1+floor((i2-i1)*r)
!
!        ipiv = i1+(i2-i1)/2
!
        ival  = ivals(ipiv)
!
!       Swap the pivot element and the last element.
!
        ivals(ipiv) = ivals(i2)
        ivals(i2)   = ival
!
        i3 = i1
!
        do 1000 i=i1,i2-1
        if( ivals(i) .lt. ival) then
        id  = ivals(i)
!
        ivals(i)  = ivals(i3)
        ivals(i3) = id
!
        i3=i3+1
        endif
 1000 continue
!
        id        = ivals(i3)
        ivals(i3) = ivals(i2)
        ivals(i2) = id
!
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
!  Reduce to upper triangular 
!
do i=1,n-1
do j=n,i+1,-1
aa(1)=a(i,j-1)
aa(2)=a(i,j)

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


function eye(n) result(a)
implicit double precision (a-h,o-z)

double precision, allocatable :: a(:,:)
integer n

!
!  Return an identity matrix which is dimensioned (n,n).
!

allocate(a(n,n))
a = 0
do i=1,n
a(i,i)=1.0d0
end do

end function

end module
