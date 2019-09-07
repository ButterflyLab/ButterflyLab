module mexprint

interface prin2
   module procedure prin2_0
   module procedure prin2_1
   module procedure prin2_2
end interface prin2

interface prini
   module procedure prini_0
   module procedure prini_1
   module procedure prini_2
end interface prini

contains

subroutine prina(str)
implicit double precision (a-h,o-z)
character(len=*)    :: str
character(len=1024) :: out
write(out,'(A,"\n")') str
call mexprintf(out)
end subroutine


subroutine prin2_0(str,a)
implicit double precision (a-h,o-z)

character(len=*), intent(in) :: str
double precision, intent(in) :: a
character(len=2048)          :: out

write (out,'(A)') str
call mexprintf(out)
call mexprintf("\n")

write (out,'(8(2x,e15.7))') a
call mexprintf(out)
call mexprintf("\n")

end subroutine


subroutine prin2_1(str,a)
implicit double precision (a-h,o-z)

character(len=*), intent(in) :: str
double precision, intent(in) :: a(:)
character(len=2048)          :: out

write (out,'(A)') str
call mexprintf(out)
call mexprintf("\n")

n   = floor(size(a,1)/5.0d0)
k   = size(a,1) - n*5
idx = 0

do i=1,n
do j=1,5
idx = idx+1
write (out,'(2x,e15.7)') a(idx)
call mexprintf(out)
end do
call mexprintf("\n")
end do

do j=1,k
idx = idx+1
write (out,'(2x,e15.7)') a(idx)
call mexprintf(out)
end do
call mexprintf("\n")

end subroutine


subroutine prin2_2(str,a0)
implicit double precision (a-h,o-z)

character(len=*), intent(in)  :: str
double precision, intent(in)  :: a0(:,:)
character(len=2048)           :: out
double precision, allocatable :: a(:)

m = size(a0,1)*size(a0,2)
allocate(a(m))
a = reshape(a0,[m])

call prin2_1(str,a)


end subroutine



subroutine prini_0(str,a)
implicit double precision (a-h,o-z)

character(len=*), intent(in) :: str
integer, intent(in)          :: a
character(len=2048)          :: out

write (out,'(A)') str
call mexprintf(out)
call mexprintf("\n")

write (out,'(2x,I8.8)') a
call mexprintf(out)
call mexprintf("\n")

end subroutine


subroutine prini_1(str,a)
implicit double precision (a-h,o-z)

character(len=*), intent(in) :: str
integer, intent(in)          :: a(:)
character(len=2048)          :: out

write (out,'(A)') str
call mexprintf(out)
call mexprintf("\n")

n   = floor(size(a,1)/5.0d0)
k   = size(a,1) - n*5
idx = 0

do i=1,n
do j=1,5
idx = idx+1
write (out,'(2x,I8.8)') a(idx)
call mexprintf(out)
end do
call mexprintf("\n")
end do

do j=1,k
idx = idx+1
write (out,'(2x,I8.8)') a(idx)
call mexprintf(out)
end do
call mexprintf("\n")

end subroutine


subroutine prini_2(str,a0)
implicit double precision (a-h,o-z)

character(len=*), intent(in)  :: str
integer, intent(in)           :: a0(:,:)
character(len=2048)           :: out
integer, allocatable          :: a(:)

m = size(a0,1)*size(a0,2)
allocate(a(m))
a = reshape(a0,[m])

call prini_1(str,a)

end subroutine

end module
