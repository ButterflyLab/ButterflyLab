module idecomp

use utils

contains

subroutine idecomp_construct(eps,a,krank,ipivs,r)
implicit double precision (a-h,o-z)

double precision               :: eps,a(:,:)
integer, allocatable           :: ipivs(:)
double precision, allocatable  :: r(:,:)

!
!
!

! double precision, allocatable :: c(:,:),q(:,:),rnorms(:),work(:),rinv(:,:),r0(:,:)
! integer, allocatable          :: idxs(:)

double precision, allocatable :: q(:,:),rnorms(:),c(:,:),rinv(:,:),work(:),r2(:,:),r0(:,:)

n   = size(a,1)
m   = size(a,2)

allocate(rnorms(n+m+1000),q(n,m),ipivs(m))

q = a
call gspiv(q,n,m,eps,rnorms,ipivs,krank)
! call prini("krank = ",krank)
! call prin2("rnorms = ",rnorms(1:krank))
! call prini("ipivs = ",ipivs)



call insorti(krank,ipivs(1:krank))
! call prini("ipivs = ",ipivs(1:krank))

allocate(r2(krank,m),c(n,m),rinv(krank,krank),work(krank*krank*8 + 1000))

c = a(:,ipivs)

r2 = matmul(transpose(q(:,1:krank)),c)
rinv = r2(1:krank,1:krank)
call orthom(rinv,krank,work,cond)

allocate(r(krank,m),r0(krank,m))
r0 = matmul(rinv,r2)

! print *,norm2(c - matmul(c(:,1:krank),r0))
! print *,norm2(a - matmul(c(:,1:krank),r(:,ipivs)))

do j=1,m
r(:,ipivs(j)) = r0(:,j)
end do

! r = 0
! do j=1,krank
! r(j,j) = 1
! end do
! r(:,krank+1:m) = matmul(rinv,c(:,krank+1:m))


! ! reorder r
! ! r = r(:,ipivs)
! do j=1,m
! r(:,j) = r0(:,ipivs(j))
! end do

end subroutine




end module
