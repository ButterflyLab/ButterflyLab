        subroutine corthom(a,n,work,cond)
        implicit double precision (a-h,o-z)
        save
        double complex a(1),work(1)
c 
c        this subroutine inverts the user-supplied matrix
c        a by means of the gram-schmidt process. in fact, this
c        is a memory management routine for the routine cortho1,
c        which performs the actual inversion
c 
c               input parameters:
c 
c  a - the n*n - matrix to be inverted (destroyed by the
c         subroutine)
c  n - the dimensionality of the matrix a
c 
c               output parameters:
c 
c  a - the inverse of a on input
c  cond - the estimate of the condition number of the matrix a
c         that has been inverted
c 
c               work array
c 
c  work - must be at least 2 * n * (n+1)+2 double precision locations long
c 
c 
c        . . . allocate memory for the matrix inverter
c 
        ib=1
        lb=n*n
c 
        iw=ib+lb+1
c 
c       invert the matrix
c 
        call cortho1(a,work(ib),n,work(iw),cond)
        return
        end
c 
c 
c 
c 
c 
        subroutine cortho1(a,b,n,work,cond)
        implicit double precision (a-h,o-z)
        save
        double complex a(n,n),b(n,n),cd,zero,
     1      work(1)
c 
c       set the matrix b to unity
c 
        done=1
        zero=0
        do 1400 i=1,n
        do 1200 j=1,n
        b(j,i)=zero
 1200 continue
        b(i,i)=1
 1400 continue
c 
c        conduct the gram-schmidt process
c 
        dnormax=-1
        dnormin=1.0d23
c 
        do 4000 i=1,n
c 
c       orthogonalize the i-th column to all preceeding ones
c 
        if(i .eq. 1) goto 2190
        do 2180 j=1,i-1
c 
        call corthom_scap(a,n,i,j,cd)
c 
        do 2170 k=1,n
        a(i,k)=a(i,k)-a(j,k)*cd
        b(i,k)=b(i,k)-b(j,k)*cd
 2170 continue
 2180 continue
 2190 continue
c 
c       normalize the i-th row
c 
        call corthom_scap(a,n,i,i,cd)
c 
        d=cd
        d=done/sqrt(d)
        if(dnormax .lt. d) dnormax=d
        if(dnormin .gt. d) dnormin=d
        do 2400 j=1,n
        a(i,j)=a(i,j)*d
        b(i,j)=b(i,j)*d
 2400 continue
c 
c       orthogonalize all subsequent rows to the i-th one
c 
        if(i .eq. n) goto 4000
c 
        do 3000 j=i+1,n
c 
        cd=0
        do 2600 k=1,n
        cd=cd+a(i,k)*dconjg(a(j,k))
 2600 continue
        cd=dconjg(cd)
c 
        do 2800 k=1,n
        a(j,k)=a(j,k)-cd*a(i,k)
        b(j,k)=b(j,k)-cd*b(i,k)
 2800 continue
 3000 continue
c 
 4000 continue
        cond=dnormax/dnormin
c 
c       now, multiply the abjoint of the resulting
c       orthogonal matrix by the triangular one,
c       obtaining the inverse of the original a
c 
        do 5000 i=1,n
        do 4800 j=1,n
        cd=0
        do 4200 k=1,n
        cd=cd+dconjg(a(k,i))*b(k,j)
 4200 continue
        work(j)=cd
 4800 continue
        do 4900 j=1,n
        a(j,i)=work(j)
 4900 continue
 5000 continue
c 
c        now, transpose a
c 
        do 5400 i=2,n
        do 5200 j=1,i
        cd=a(i,j)
        a(i,j)=a(j,i)
        a(j,i)=cd
 5200 continue
 5400 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine corthom_matvec(a,n,x,y)
        implicit double precision (a-h,o-z)
        save
        double complex a(n,n),x(1),y(1),cd
c 
c        apply the matrix a to the vector x obtaining y
c 
        do 1400 i=1,n
        cd=0
        do 1200 j=1,n
        cd=cd+a(i,j)*x(j)
 1200 continue
        y(i)=cd
 1400 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine corthom_scap(a,n,i,j,prod)
        implicit double complex (a-h,o-z)
        save
        dimension a(n,n)
c 
        prod=0
        do 1200 k=1,n
        prod=prod+a(i,k)*conjg(a(j,k))
 1200 continue
        return
        end
  
  
