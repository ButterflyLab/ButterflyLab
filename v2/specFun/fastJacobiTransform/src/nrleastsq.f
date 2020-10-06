        implicit real *8 (a-h,o-z)
        dimension a(50 000),c(50 000),c2(50 000),
     1      diff(50 000),w(1000 000),proj(1000 000)
c 
        call prini(6,13)
C 
C       SET ALL PARAMETERS
C 
         PRINT *, 'ENTER n'
         READ *,n
         CALL PRINF('n=*',n,1 )
C 
         PRINT *, 'ENTER m'
         READ *,m
         CALL PRINF('m=*',m,1 )
c 
C 
         PRINT *, 'ENTER k'
         READ *,k
         CALL PRINF('k=*',k,1 )
c 
c 
        PRINT *, 'ENTER -log(eps) base 10'
        READ *,db
        CALL PRIN2('db=*',db,1 )
c 
        eps=0.1**db
        call prin2('and eps=*',eps,1)
c 
c       test the gram-schmidt procedure
c 
        lw=300000
        call leastes(a,n,m,k,eps,c,c2,diff)  
        
        call nrleas_proj(a,m,k,eps,proj,ncols,w)
  
  
  
  
        call nrleas_proj_test(a,n,m,eps,proj,w)

  
        stop
        end
c 
c 
c 
c 
c 
        subroutine nrleas_proj_test(a,n,m,eps,proj,w)
        implicit real *8 (a-h,o-z)
        save
        real *8 a(n,m),w(n,m),proj(n,m),
     1      uu(1000 000),vv(1000 000),ss(1000)
c
        call prin2('and proj=*',proj,n*n)
c
c        check that proj is in fact a projection
c
        call rleamat(ier,proj,n,n,a,n,n,w)

        call prin2('and w=*',w,n*n)
c
        do 3400 i=1,n
        do 3200 j=1,n
c
        w(j,i)=w(j,i)-a(j,i)
 3200 continue
 3400 continue

        call prin2('and differences*',w,n*n)
c
c       finally, look at singular values
c
        eps2=1.0d-16
        lw=1000 000
        call svdpivot(ier,proj,n,n,uu,vv,ss,ncols2,eps2,
     1      w,lw,ltot)

        call prin2('and ss=*',ss,ncols2)
        call prinf('and ncols2=*',ncols2,1)

        return
        end
c 
c 
c 
c 
c 
        subroutine leastes(a,n,m,k,eps,c,c2,diff)
        implicit real *8 (a-h,o-z)
        save
        real *8 a(m,k),w(500 000),
     1      c(n,k),rnorms(1000),tt(1000),
     2      c2(n,k),ima,x(100000),diff(n,k)
c 
        data ima/(0.0d0,1.0d0)/
c 
c        construct the test matrices a, c
c 
        h=1
        h=h/n
        do 1100 i=1,2000
        tt(i)=i*h+ima/10
cccc        tt(i)=i*h
 1100 continue
c 
        do 1400 i=1,k
        do 1200 j=1,m
        a(j,i)=exp(-tt(i)*tt(j))+exp(-tt(i)/tt(j))
c 
cccc        if(i .eq. j) a(j,i)=a(j,i)+10
c 
 1200 continue
 1400 continue
  
  
  
        call prin2('a as created*',a,n*m)
  
c 
        do 1800 i=1,k
        do 1600 j=1,n
        c(j,i)=exp(-tt(i)*tt(j))+exp(-tt(i)/tt(j))
c 
cccc        if(i .eq. j) c(j,i)=c(j,i)+10
  
 1600 continue
 1800 continue
c 
cccc        call prin2('c as created*',c,n*k)
  
c 
c        solve the equation A * X = C in the least squares sense
c 
        ifcheck=1
  
  
        t1=clotatim()
        call nrleamatrr(a,n,m,k,c,x,eps,ncols,
     1    rnorms,w,ifcheck,errl2,errmax)
  
  
        t2=clotatim()
  
        call prin2('and time=*',t2-t1,1)
c 
cccc        call prin2('x as calculated*',x,n*k)
        call prin2('after nrleamatrr, errl2=*',errl2,1)
        call prin2('after nrleamatrr, errmax=*',errmax,1)
cccc        call prin2('after nrleamatrr, w=*',w,n*k)
  
  
  
       return
       end
c 
c 
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c       this is the end of the debugging code and the beginning of the
c       least squares subroutines proper
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c         This file contains 8 user-callable subroutines: nrleastsq,
c         nrleastsq2, nrleamatlr, rleamatr, nrleamatrr, nrleamatll,
c         nrleamatrl, nrleas_proj. Following is a brief description 
c         of these subroutines.
c
c   nrleas_proj - constructs the orthogonal projection operator from 
c         R^n to R^n, projecting R^n on its subspace spanned by the 
c         collection of user-specified vectors
c 
c    NOTE: PLEASE NOTE THAT THE SUBROUTINE nrleamatrl IS ALMOST NOT
c         A USER-CALLABLE ONE. WHILE IT IS CONCIEVABLE THAT SOME
c         USERS WILL FIND IT USEFUL, MOST ARE EXPECTED TO BE HAPPY
c         WITH THE SUBROUTINES nrleamatrr, nrleamatll, TO THE EXTENT
c         THAT THEY ARE HAPPY WITH ANYTHING IN THIS FILE.
c 
c   nrleastsq - constructs a QR-type (not exactly QR) decomposition
c         of the real input matrix A, to be used by the subroutine
c         nrleastsq2 for the solution of least squares problems of
c         the form
c 
c               A X \sim Y,                                             (1)
c 
c        and by the subroutine nrleamatrl for the least squares
c        solution of the matrix equation
c 
c         A(n,m) * X(m,k) = C(n,k).                                     (2)
c 
c   nrleastsq2 - uses the QR-type (not exactly QR) decomposition
c         of the matrix A produced by a prior call to nrleastsq (see)
c         to solve in the least squares sense the linear system
c 
c                    A X = RHS.                                         (3)
c 
c   rleamatr - solves in the least squares sense the matrix
c        equation
c 
c                     A(k,l) * X(l,m) * B(m,n) = C(k,n),                (3a)
c 
c       where A, B, C are user-specified real matrices, and X
c       is the matrix to be found. Note that the dimensionalities
c       of the matrices A,B,C, X are as general as could be
c 
c   nrleamatlr - solves in the least squares sense the matrix
c        equation
c 
c                     A(k,l) * X(l,m) * B(m,n) = C(k,n),                (3b)
c 
c       where A, B, C are user-specified real matrices, and X
c       is the matrix to be found. Note that the dimensionalities
c       of the matrices A,B,C, X are as general as could be
c 
c       PLEASE NOTE THAT THIS SUBROUTINE CONSTRUCTS GENERALIZED
c       INVERSES OF THE MATRICES A, B INDEPENDENTLY FROM EACH OTHER.
C       THUS, A CONSPIRACY IS POSSIBLE THAT WILL CAUSE SOME OF THE
C       ELEMENTS OF THE MATRIX X TO BE OF THE ORDER ||C||/EPS^2,
C       AS OPPOSED TO THE PROPER ||C||/EPS. CLEARLY, THIS WILL NOT
C       HAPPEN IF THE MATRICES A,B (WHILE RECTANGULAR) HAVE NO SMALL
C       SINGULAR VALUES. THE PURPOSE OF THIS SHORT-CUT IS TO AVOID
C       THE CONSTRUCTION OF SINGULAR VALUE DECOMPOSITIONS OF A, B;
C       THE LATTER PROCEDURE IS KNOWN TO BE POTENTIALLY EXPENSIVE
C       AND EVEN UNRELIABLE. THE PROPER PROCEDURE USING THE SVDs
C       IS PERFORMED BY THE SUBROUTINE RLEAMATR (SEE)
c 
c   nrlearmatll - solves in the least squares sense the matrix
c         equation
c 
c                    A(n,m) * X(m,k) = C(n,k),                          (4)
c 
c         without using any additional data (i.e. it performs all
c         factorizations itself).
c 
c   nrlearmatrr - solves in the least squares sense the matrix
c         equation
c 
c                    X(n,m) * A(m,k) = C(n,k),                          (5)
c 
c         without using any additional data (i.e. it performs all
c         factorizations itself).
c 
c   nrleamatrl - solves in the least squares sense the
c        matrix equation
c 
c                    A(n,m) * X(m,k) = C(n,k),                          (6)
c 
c        using as input the matrix c and the array w obtained
c        via a preceding call to the subroutine nrleastsq
c        (see).
c 
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c 
c 
c 
c 
c 
c 
        subroutine nrleas_proj(a,n,m,eps,proj,ncols,w)
        implicit real *8 (a-h,o-z)
        save
        real *8 a(n,m),w(n,m),cd,proj(n,m)
c
c       This subroutine constructs the orthogonal projection 
c       operator from R^n to R^n, projecting R^n on its subspace
c       spanned by the collection of user-specified vectors
c
c                  Input parameters:
c
c  a - collection of m real vectors of length n
c  n - the length of the vectors (dimensionality of the space)
c  m - the number of vectors (dimensionality of subspace)
c  eps - the precision of calculation
c  
c                  Output parameters:
c
c  proj - the matrix of orthogonal projection on the on the 
c       subspace in R^n spanned by the columns of the matrix a
c  ncols - the rank of the matrix proj
c
c                  work arrays:
c
c  w - must be at least n*(n+1) real *8 elements long
c 
c        . . . gram-schmidt the columns of a
c
        call rleamem(a,w,n*m)
c
        ifpivot=1
        call rleasgrm(w,n,m,w(1,m+1),eps,ncols,ifpivot)
c
c        now, construct the projection matrix
c
        do 2600 i=1,n
        do 2400 j=1,n
c
        cd=0
        do 2200 ll=1,ncols
c
        cd=cd+w(i,ll)*w(j,ll)
 2200 continue
c
        proj(i,j)=cd
 2400 continue
 2600 continue
c
        return
        end
c
c
c
c
c
        subroutine nrleamatlr(ier,a,b,c,k,l,m,n,eps,x,
     1      w,lw,ltot)
        implicit real *8 (a-h,o-z)
        save
        real *8  a(k,l),x(l,m),b(m,n),c(k,n),w(1)
c 
c        This subroutine solves in the least squares sense the
c        matrix equation
c 
c                     A(k,l) * X(l,m) * B(m,n) = C(k,n),               (1)
c 
c       where A, B, C are user-specified real matrices, and X
c       is the matrix to be found. Note that the dimensionalities
c       of the matrices A,B,C, X are as general as could be
c 
c       PLEASE NOTE THAT THIS SUBROUTINE CONSTRUCTS GENERALIZED
c       INVERSES OF THE MATRICES A, B INDEPENDENTLY FROM EACH OTHER.
C       THUS, A CONSPIRACY IS POSSIBLE THAT WILL CAUSE SOME OF THE
C       ELEMENTS OF THE MATRIX X TO BE OF THE ORDER ||C||/EPS^2,
C       AS OPPOSED TO THE PROPER ||C||/EPS. CLEARLY, THIS WILL NOT
C       HAPPEN IF THE MATRICES A,B (WHILE RECTANGULAR) HAVE NO SMALL
C       SINGULAR VALUES. THE PURPOSE OF THIS SHORT-CUT IS TO AVOID
C       THE CONSTRUCTION OF SINGULAR VALUE DECOMPOSITIONS OF A, B;
C       THE LATTER PROCEDURE IS KNOWN TO BE POTENTIALLY EXPENSIVE
C       AND EVEN UNRELIABLE. THE PROPER PROCEDURE USING THE SVDs
C       IS PERFORMED BY THE SUBROUTINE RLEAMATR (SEE)
c 
c                          input parameters:
c 
c  a,b,c - matrices in (1)
c  k,l,m,n - dimensionalities in (1)
c  eps - the parameter telling the subroutine at which point a
c       a singular value of the matrix is to be declared to be zero
c       and ignored (see subroutine rleasas for details)
c  lw - the length of the user-supplied work array
c 
c                          output parameters:
c 
c  ier - error return code
c      ier=0 means successful execution
c      ier .ne. 0 means that the amount of space lw allocated supplied
c            in the work array w is insufficient
c  x - the unknown matrix in (1)
c  ltot - the length of the work array w actually used by the subroutine
c 
c                          work array:
c 
c  w - must be sufficiently long
c 
c 
c        . . . start with finding the matrix (approximately)
c              satisfying the equation
c 
c                     A(k,l) * Y(l,n) = C(k,n)                         (2)
c 
        ier=0
c 
        iy=1
        ly=l*n+100
c 
        irnorms=iy+ly
        lr=k
        if(lr .lt. l) lr=l
        if(lr .lt. m) lr=m
        if(lr .lt. n) lr=n
c 
        lrnorms=lr
c 
        iw=irnorms+lrnorms
        lenw=n*m
        if(n*k .gt. lenw) lenw=n*k
        if(m*k .gt. lenw) lenw=m*k
c 
        len2=lenw*8+500
c 
        ltot=iw+lenw
        if(ltot .gt. lw) then
            ier=16
            return
        endif
c 
        ifcheck=0
        call nrleamatll(a,k,l,n,c,w(iy),eps,ncols,
     1    w(irnorms),w(iw),ifcheck,errl2,errmax)
c 
cccc        call prin2('in rleamatlr after rleamatll, w(irnorms)=*',
cccc     1      w(irnorms),ncols)
c 
c        now, solve the equation
c 
c        X(l,m) * B(m,n) = Y(l,n)                                      (3)
c 
        call nrleamatrr(b,l,m,n,w(iy),x,eps,ncols2,
     1    w(irnorms),w(iw),ifcheck,errl2,errmax)
  
cccc        call prin2('in rleamatlr after rleamatrr, w(irnorms)=*',
cccc     1      w(irnorms),ncols2)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine rleamatr(ier,a,b,c,k,l,m,n,eps,delta,x,
     1      w,lw,ltot)
        implicit real *8 (a-h,o-z)
c 
        save
        real *8  a(k,l),x(l,m),b(m,n),c(k,n),w(1)
c 
c        this subroutine solves in the least squares sense the
c        equation
c 
c                     A(k,l) * X(l,m) * B(m,n) = C(k,n),               (1)
c 
c       where A, B, C are user-specified real matrices, and X
c       is the matrix to be found. Note that the dimensionalities
c       of the matrices A,B,C, X are as general as could be
c 
c 
c                          input parameters:
c 
c  a,b,c - matrices in (1)
c  k,l,m,n - dimensionalities in (1)
c  eps - the parameter telling the subroutine at which point a
c       a singular value of the matrix is to be declared to be zero
c       and ignored (see subroutine rleasas for details)
c  delta - the "Tikhonov constant" (see subroutine rleasas for details)
c  lw - the length of the user-supplied work array
c 
c                          output parameters:
c 
c  ier - error return code
c      ier=0 means successful execution
c      ier .ne. 0 means that the amount of space lw allocated supplied
c            in the work array w is insufficient
c  x - the unknown matrix in (1)
c  ltot - the length of the work array w actually used by the subroutine
c 
c                          work array:
c 
c  w - must be sufficiently long
c 
c 
c        . . . construct the SVD of A
c 
        ier=0
        isa=1
        lsa=min(k,l)+2
c 
        iua=isa+lsa
        lua=k*l+10
c 
        iva=iua+lua
        lva=k*l+10
c 
        iwork=iva+lva
c 
        eps2=1.0d-13
        lw1=lw-iwork-2
c 
        call prinf('in rleamatr, lw1=*',lw1,1)
c 
  
cccc        subroutine svdpivot(ier,a,n,m,u,v,s,ncols,eps,
cccc     1      w,lw,ltot)
  
  
        call svdpivot(ier,a,k,l,w(iua),w(iva),w(isa),ncolsa,eps2,
     1      w(iwork),lw1,ltot1)
c 
        ltot=ltot1+iwork
c 
        if(ier .ne. 0) return
c 
        call prinf('after first svdpivot, ncolsa=*',ncolsa,1)
        call prin2('after first svdpivot, sa=*',w(isa),ncolsa)
c 
c        perform garbage collection
c 
        lua2=k*ncolsa+2
        lva2=l*ncolsa+2
        iva2=iua+lua2
        lua=lua2
c 
        call rleamem(w(iva),w(iva2),lva2)
c 
        iva=iva2
        lva=lva2
c 
c        . . . construct the SVD of B
c 
        isb=iva+lva
        lsb=min(m,n)+2
c 
        iub=isb+lsb
        lub=m*n+10
c 
        ivb=iub+lub
        lvb=m*n+10
c 
        iwork=ivb+lvb
        lw2=lw-iwork
cccc        call prinf('in rleamatr, lw2=*',lw2,1)
c 
        eps2=1.0d-13
c 
        call svdpivot(ier,b,m,n,w(iub),w(ivb),w(isb),ncolsb,eps2,
     1      w(iwork),lw2,ltot2)
c 
        ii=iwork+ltot2
        if(ii .gt. ltot) ltot=ii
c 
cccc        call prinf('after second svdpivot, ltot=*',ltot,1)
c 
        if(ier .ne. 0) return
c 
        call prinf('after second svdpivot, ncolsb=*',ncolsb,1)
        call prin2('after second svdpivot, sb=*',w(isb),ncolsb)
c 
c        perform garbage collection
c 
        lub2=m*ncolsb+2
        lvb2=n*ncolsb+2
        ivb2=iub+lub2
c 
        call rleamem(w(ivb),w(ivb2),lvb2)
c 
        ivb=ivb2
        lvb=lvb2
c 
        iwork=ivb+lvb
        lwork=ncolsa*n+10
c 
        ii=iwork+lwork
        call prinf('in rleamatr, final ii=*',ii,1)
        if(ii .gt. ltot) ltot=ii
        if(ltot .lt. lw) goto 2200
        ier=2
        return
 2200 continue
c 
c         perform the remainder of the solution
c 
        call rleamat0(c,k,l,m,n,eps,delta,x,
     1      w(isa),w(iua),w(iva),w(isb),w(iub),w(ivb),
     2      w(iwork),lw,ltot,ncolsa,ncolsb)
  
  
        return
        end
  
c 
c 
c 
c 
c 
        subroutine nrleamatrr(a,n,m,k,c,x,eps,ncols,
     1    rnorms,w,ifcheck,errl2,errmax)
        implicit real *8 (a-h,o-z)
        save
        dimension w(1)
        real *8 x(n,m),c(n,k),a(m,k)
c 
c        This subroutine solves in the least squares sense the
c        matrix equation
c 
c        X(n,m) * A(m,k) = C(n,k).                              (1)
c 
c                     input parameters:
c 
c  a - the matrix to be decomposed (note that it is NOT damaged
c       by this subroutine in any way)
c  n, m, k - the dimensions in (1) above
c  c - the right-hand side above NOT damaged by this subroutine
c  eps - the accuracy to which the decomposition (1) will be computed.
c  ifcheck - integer parameter telling the subroutine whether the
c      error of the approximations should be calculated.
c   ifcheck=1 will cause the subroutine to evaluate (and return)
c      the errors - both maximum and relative l^2
c   ifcheck=0 will cause the subroutine skip the error evaluation.
c      This is a CPU time saving feature.
c 
c                     output parameters:
c 
c  x - the solution of (1) above in the least squares sense
c  ncols - the rank of the matrix a to the precision eps. also the
c        second dimension of the matrices u, w on exit
c  rnorms - the normalizing factors in the gram-schmidt process.
c        only the first ncols of them are meaninful, but the array
c        has to be dimensioned by the user to be at least m+1
c        real *8 elements long.
c  w - the first n*k real *8 elements of w contain the
c        discrepancies in the approximations of the n*k elements
c        of the matrix c
c  errl2 - the l^2 error of the obtained approximation; in other
c        words, errl2=sqrt(sum_{i=1}^{n*k} |w_i|^2)
c 
c                    work arrays:
c 
c  w - must be at least n*m*16 + 500 real *8 elements long
c 
c 
c       . . . transpose both the matrices a, c
c 
        call rleastra(a,m,k,w)
        call rleascop(w,a,m*k)
c 
        call rleastra(c,n,k,w)
        call rleascop(w,c,n*k)
c 
c       solve least-squares problem
c 
c       A^* (k,m) * x^* (m,n) = C^* (k,n)
c 
        call nrleamatll(a,k,m,n,c,x,eps,ncols,
     1    rnorms,w,ifcheck,errl2,errmax)
c 
c       transpose everything back
c 
        iaia=n*m
        if(k*m .gt. iaia) iaia=k*m
        if(k*n .gt. iaia) iaia=k*n
  
        iaia=iaia*2+20
  
        call rleascop(w,w(iaia),n*k)
c 
        call rleastra(a,k,m,w)
        call rleascop(w,a,m*k)
c 
        call rleastra(c,k,n,w)
        call rleascop(w,c,n*k)
c 
        call rleastra(x,m,n,w)
        call rleascop(w,x,n*m)
c 
        call rleastra(w(iaia),k,n,w)
  
        return
        end
c 
c 
c 
c 
        subroutine nrleamatll(a,n,m,k,c,x,eps,ncols,
     1    rnorms,w,ifcheck,errl2,errmax)
        implicit real *8 (a-h,o-z)
        save
        real *8 x(m,k),c(n,k),a(n,m),cd,w(1),cd2
c 
c        This subroutine solves in the least squares sense the
c        matrix equation
c 
c         A(n,m) * X(m,k) = C(n,k).                            (1)
c 
c                     input parameters:
c 
c  a - the matrix to be decomposed (note that it is NOT damaged
c       by this subroutine in any way)
c  n, m, k - the dimensions in (1) above
c  c - the right-hand side above, NOT damaged by this subroutine
c  eps - the accuracy to which the decomposition (1) will be computed.
c  ifcheck - integer parameter telling the subroutine whether the
c       error of the approximations should be calculated.
c    ifcheck=1 will cause the subroutine to evaluate (and return)
c       the errors - both maximum and relative l^2
c    ifcheck=0 will cause the subroutine skip the error evaluation.
c       This is a CPU time saving feature.
c 
c                     output parameters:
c 
c  x - the solution of (1) above in the least squares sense
c  ncols - the rank of the matrix a to the precision eps. also the
c        second dimension of the matrices u, w on exit
c  rnorms - the normalizing factors in the gram-schmidt process.
c        only the first ncols of them are meaninful, but the array
c        has to be dimensioned by the user to be at least m+1
c        real *8 elements long.
c  w - the first n*k real *8 elements of w contain the
c        discrepancies in the approximations of the n*k elements
c        of the matrix c
c  errl2 - the l^2 error of the obtained approximation; in other
c        words, errl2=sqrt(sum_{i=1}^{n*k} |w_i|^2)
c 
c                    work arrays:
c 
c  w - must be at least n*m*8 + 500 real *8 elements long
c 
c 
c        . . . decompose the matrix a
c 
        call nrleastsq(a,n,m,eps,ncols,rnorms,w)
c 
c       obtain the matrix x
c 
        call nrleamatrl(w,c,x,n,m,k)
c 
c        multiply A by X and check how close the result is to C
c 
        if(ifcheck .eq. 0) return
  
        call nrleamult(a,x,w,n,m,k)
        call rleasubt(w,c,w,n*k)
c 
        call rleascap(w,w,n*k,cd)
        call rleascap(c,c,n*k,cd2)
c 
        errl2=dsqrt(abs(cd/cd2))
  
        errmax=0
        do 2200 i=1,n*k
c 
        d=w(i)*(w(i))
        if(errmax .lt. d) errmax=d
 2200 continue
c 
        errmax=dsqrt(errmax)
c 
        return
        end
  
  
  
c 
c 
c 
c 
c 
        subroutine nrleamatrl(w,c,x,n,m,k)
        implicit real *8 (a-h,o-z)
        save
        dimension w(1)
        real *8 x(m,k),c(n,k)
c 
c        This subroutine solves in the least squares sense the
c        matrix equation
c 
c         A(n,m) * X(m,k) = C(n,k),                            (1)
c 
c        using as input the matrix c and the array w obtained
c        via a preceding call to the subroutine nrleastsq
c        (see).
c 
c                     input parameters:
c 
c  w - the array containing the decomposition of the matrix a,
c        obtained via a preceding call to nrleastsq  (see)
c  c - the right-hand side above
c  n, m, k - the dimensions in (1) above
c 
c                     output parameters:
c 
c  x - the solution of (1) above in the least squares sense
c 
c 
c        . . . solve the least squares problems for columns of x
c              one after another
c 
        do 1400 i=1,k
c 
        call nrleasts2(w,c(1,i),x(1,i) )
 1400 continue
  
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine nrleastsq(a,n,m,eps,ncols,rnorms,w)
        implicit real *8 (a-h,o-z)
        save
        real *8 a(n,m)
        dimension rnorms(1),w(1)
c 
c        This subroutine constructs the decomposition of the input
c        matrix a, to be used by the subroutine nrleastsq2 for the
c        solution of least squares problems of the form
c 
c               A X \sim Y,                                          (1)
c 
c        and by the subroutine nrleamatrl for the least squares
c        solution of the matrix equation
c 
c 
c         A(n,m) * X(m,k) = C(n,k).                                  (2)
c 
c        The decomposition is stored in the array w; the first
c        n*m*8 + 500 real *8 elements of the latter should not be
c        changed between the call to this subroutine and the
c        subsequent calls to nrleastsq2.
c 
c   NOTE: this subroutine uses the subroutine rleastsq to perform
c        all of the work; this is simply a memory manager for
c        rleastsq.
c 
c                     input parameters:
c 
c  a - the matrix to be decomposed (note that a is NOT damaged
c       by this subroutine in any way)
c  n, m - the dimensions of a
c  eps - the accuracy to which the decomposition (1) will be computed.
c 
c                     output parameters:
c 
c  w - the array containing the decomposition to be used by nrleastsq2.
c         Must be at least n*m*8 + 500 real *8 elements long
c  rnorms - the normalizing factors in the gram-schmidt process.
c        only the first ncols of them are meaninful, but the array
c        has to be dimensioned by the user to be at least m+1
c        real *8 elements long.
c 
c 
c        . . . allocate memory for the decomposition of the matrix
c 
        iu=20
        lu=(n*m+10)*2
c 
        iw=iu+lu
        lw=(n*m+10)*2
c 
        it=iw+lw
        lt=(n*m+10)*2
c 
        iv=it+lt
        lv=(n*m+10)*2
c 
c       decompose the matrix
c 
        call rleastsq(a,w(iu),w(iw),w(it),n,m,ncols,
     1      rnorms,eps,w(iv))
c 
c       store various data at the beginning of the array w
c 
        w(1)=n+0.1
        w(2)=m+0.1
        w(3)=ncols+0.1
        w(4)=eps
        w(5)=iu+0.1
        w(6)=iw+0.1
        w(7)=it+0.1
        w(8)=iv+0.1
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine nrleasts2(w,rhs,x)
        implicit real *8 (a-h,o-z)
        save
        dimension w(1)
        real *8 x(1),rhs(1)
c 
c         This subroutine uses the QR-type (not exactly QR) decomposition
c         of a matrix to solve in the least squares sense the linear
c         system
c 
c                    A X = RHS                                             (1)
c 
c         The expansion used must have been produced by a prior call
c         to the subroutine nrleastsq (see), and is supplied to
c         this subroutine via the input array w.
c 
c                     input parameters:
c 
c  w - the array produced via a preceding call to the subroutine
c         nrleastsq (see); please note that the first  n*m*8 + 500
c         real *8 elements of this array should not be changed
c         between the call to this subroutine and the preceding call
c         to nrleastsq
c  rhs - the right-hand side in (1) (of length n)
c 
c                     output parameters:
c 
c  x - the solution (of length m) of the system (1) in the
c         least squares sense
c 
c 
c        . . . retrieve from the beginning of the array w all of the
c              relevant integer data
c 
        n=w(1)
        m=w(2)
        ncols=w(3)
        iu=w(5)
        iw=w(6)
        it=w(7)
        iv=w(8)
c 
c       solve the least squares problem
c 
        call rleasts2(w(iu),w(iw),w(it),n,m,ncols,rhs,
     1       x,w(iv) )
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine rleastsq(a,u,w,t,n,m,ncols,rnorms,eps,v)
        implicit real *8 (a-h,o-z)
        save
        real *8 a(n,m),u(n,1),v(m,1),w(m,1),t(1)
        dimension rnorms(1)
c 
c        this subroutine decomposes the input matrix a in the
c        form
c               A=U * T * W^*                                       (1)
c 
c        with the columns of each of the matrices U,W orthonormal,
c        and T a lower triangular matrix of dimensionality (ncols,ncols)
c 
c   NOTE: this subroutine uses the subroutine leastsq1 (see) to perform
c        almost all work. However, when m < n, it tends to be much
c        more efficient than leastsq1, since it performs the two
c        Gram-Schmidt processes in the optimal order (starting with
c        the longer gram-schmidt involving shorter vectors)
c 
c 
c                     input parameters:
c 
c  a - the matrix to be decomposed (note that a is NOT damaged
c       by this subroutine in any way)
c  n,m - the dimensionalities of a
c  eps - the accuracy to which the decomposition (1) will be computed.
c 
c                     output parameters:
c 
c  u - the matrix in (1). note that the matrix u has to be dimensioned
c        u(n,m), but on exit from this subroutine only the first
c        ncols columns of  u  are meaningful, the rest making no
c        sense whatsoever, so that effectively on exit u is dimensioned
c        u(n,ncols)
c  w - the third matrix in (1). it is dimensioned w(m,ncols) on
c        exit. on entry the parameter ncols is unknown,
c        so w should be dimensioned w(m,n) by the user. however,
c        on exit only the first ncols columns of w are meaningful
c  t - the second matrix in (1); On exit, it is structured as an
c        ncols * ncols matrix; however, on entry ncols is not
c        known, so it should be at least m*n real *8 locations long.
c  ncols - the rank of the matrix a to the precision eps. also the
c        second dimension of the matrices u, w on exit
c  rnorms - the normalizing factors in the gram-schmidt process.
c        only the first ncols of them are meaninful, but the array
c        has to be dimensioned by the user to be at least m+1
c        real *8 elements long.
c 
c 
c                      work arrays:
c 
c  v - must be at least n*m real *8 locations long
c 
c       . . . if m > n , construct the decomposition of the matrix as
c             specified by the user
c 
        if(m .lt. n-2) goto 2000
c 
        call rleast1(a,u,w,t,n,m,ncols,rnorms,eps,v)
cccc        call prin2('after first leastsq1, u=*',u,n*ncols)
c 
        return
c 
 2000 continue
c 
c       n is greater than m. transpose the matrix, and decompose
c       the transpose
c 
        call rleastra(a,n,m,v)
        call rleascop(v,a,n*m)
c 
        call rleast1(a,w,u,t,m,n,ncols,rnorms,eps,v)
c 
cccc        call prin2('after leastsq1, u=*',u,n*ncols)
c 
c        transpose back everything that needs to be transposed back
c 
        call rleastra(a,m,n,v)
        call rleascop(v,a,n*m)
c 
        call rleastra(t,ncols,ncols,v)
        call rleascop(v,t,ncols**2)
c 
        call rleasrer(t,ncols,ncols,v)
        call rleascop(v,t,ncols**2)
c 
        call rleasrec(u,n,ncols,v)
        call rleascop(v,u,n*ncols)
c 
c 
c 
        call rleasrec(w,m,ncols,v)
        call rleascop(v,w,m*ncols)
c 
        call rleasrec(t,ncols,ncols,v)
        call rleascop(v,t,ncols**2)
c 
        return
        end
  
c 
c 
c 
c 
c 
        subroutine rleasts2(u,w,t,n,m,ncols,y,x,work)
        implicit real *8 (a-h,o-z)
        save
        real *8 u(n,ncols),w(m,ncols),t(ncols,ncols)
        real *8 x(m),y(n),work(ncols),d
c 
c         this subroutine uses the QR-type (not exactly QR) decomposition
c         of a matrix to solve in the least squares sense the linear
c         system
c 
c                    A X = Y                                             (1)
c 
c         The expansion used must have been produced by a prior call
c         to the subroutine leastsq1 (see), and is of the form
c 
c               A=U * T * W^*                                            (2)
c 
c        with the columns of each of the matrices U,W orthonormal,
c        and T a lower triangular matrix of dimensionality (ncols,ncols)
c 
c                     input parameters:
c 
c  u - the matrix in (1). note that the matrix u has to be dimensioned
c        u(n,m), but on exit from this subroutine only the first
c        ncols columns of  u  are meaningful, the rest making no
c        sense whatsoever, so that effectively on exit u is dimensioned
c        u(n,ncols)
c  w - the third matrix in (1). it is dimensioned w(m,ncols) on
c        exit. on entry the parameter ncols is unknown,
c        so w should be dimensioned w(m,n) by the user. however,
c        on exit only the first ncols columns of w are meaningful
c  t - the second matrix in (1)
c  n,m - the dimensionalities of a
c  ncols - the rank of the matrix a to the precision eps. also the
c        second dimension of the matrices u, w on exit
c  y - the right-hand side in (1)
c 
c                     output parameters:
c 
c  x - the solution of the system (1) in the least squares sense
c 
c                      work arrays:
c 
c  work - must be at least ncols real *8 locations long
c 
c        . . . apply to the right-hand side the matrux  U^*
c 
        do 1400 i=1,ncols
        d=0
        do 1200 j=1,n
cccccc        d=d+u(j,i)*y(j)
        d=d+(u(j,i))*y(j)
 1200 continue
        work(i)=d
 1400 continue
c 
c       apply to the vector work the inverse of the matrix t
c 
        x(1)=work(1)/t(1,1)
        do 2000 i=2,ncols
c 
        d=0
        do 1600 j=1,i-1
        d=d+t(i,j)*x(j)
 1600 continue
c 
        x(i)=(work(i)-d)/t(i,i)
 2000 continue
c 
        do 2200 i=1,ncols
        work(i)=x(i)
 2200 continue
c 
c       multiply work by the matrix w
c 
        do 2600 i=1,m
        d=0
        do 2400 j=1,ncols
        d=d+w(i,j)*work(j)
 2400 continue
        x(i)=d
 2600 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine rleastra(a,n,m,b)
        implicit real *8 (a-h,o-z)
        save
        real *8 a(n,m),b(m,n),x(1),y(1),z(1)
c 
c       transpose a
c 
        do 1400 i=1,n
        do 1200 j=1,m
cccc        b(j,i)=a(i,j)
        b(j,i)=(a(i,j))
 1200 continue
 1400 continue
        return
c 
c 
c 
c 
        entry rleascop(x,y,n)
c 
        do 2400 i=1,n
        y(i)=x(i)
 2400 continue
        return
c 
c 
c 
c 
        entry rleasubt(x,y,z,n)
c 
        do 3400 i=1,n
        z(i)=x(i)-y(i)
 3400 continue
        return
c 
c 
c 
c 
        entry rleaadd(x,y,z,n)
c 
        do 4400 i=1,n
        z(i)=x(i)+y(i)
 4400 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine rleast1(a,u,w,t,n,m,ncols,rnorms,
     1    eps,v)
        implicit real *8 (a-h,o-z)
        save
        real *8 a(n,m),u(n,1),v(m,1),w(m,1),t(1),rnorms(1)
c 
c        this subroutine decomposes the input matrix a in the
c        form
c               A=U * T * W^*                                       (1)
c 
c        with the columns of each of the matrices U,W orthonormal,
c        and T a lower triangular matrix of dimensionality (ncols,ncols)
c 
c 
c                     input parameters:
c 
c  a - the matrix to be decomposed (note that a is NOT damaged
c       by this subroutine in any way)
c  n,m - the dimensionalities of a
c  eps - the accuracy to which the decomposition (1) will be computed.
c 
c                     output parameters:
c 
c  u - the matrix in (1). note that the matrix u has to be dimensioned
c        u(n,m), but on exit from this subroutine only the first
c        ncols columns of  u  are meaningful, the rest making no
c        sense whatsoever, so that effectively on exit u is dimensioned
c        u(n,ncols)
c  w - the third matrix in (1). it is dimensioned w(m,ncols) on
c        exit. on entry the parameter ncols is unknown,
c        so w should be dimensioned w(m,n) by the user. however,
c        on exit only the first ncols columns of w are meaningful
c  t - the second matrix in (1)
c  ncols - the rank of the matrix a to the precision eps. also the
c        second dimension of the matrices u, w on exit
c  rnorms - the normalizing factors in the gram-schmidt process.
c        only the first ncols of them are meaninful, but the array
c        has to be dimensioned by the user to be at least m+1
c        real *8 elements long.
c 
c 
c                      work arrays:
c 
c  v - must be at least n*m real *8 locations long
c 
c        . . . using gram-schmidt process with pivoting, decompose
c              the matrix a in the form
c 
c          a=U  V^*,
c 
c        with  u an orthogonal matrix of minimum rank,
c        and v whatever it wants to be
c 
        ifpivot=1
        call rleaspiv(a,n,m,u,v,ncols,rnorms,eps,ifpivot)
c 
c        using gram-schmidt process without pivoting, decompose
c        the matrix v in the form
c 
c          v=w t^*,
c 
c        with  w an orthogonal matrix of minimum rank,
c        and t a triangular matrix of dimensionality ncols * ncols
c 
        ifpivot=0
        call rleaspiv(v,m,ncols,w,t,ncols2,rnorms,eps,ifpivot)
c 
        return
        end
c 
c 
c 
c 
c 
c 
        subroutine rleaspiv(a,n,m,b,v,ncols,rnorms,eps,ifpivot)
        implicit real *8 (a-h,o-z)
        save
        real *8 a(n,m),b(n,m),v(m,1),prod
        dimension rnorms(1)
c 
c       this matrix applies the compressing gram-schmidt
c       process to the matrix a, obtaining its decomposition
c       in the form
c 
c             a=b v^T,                                        (1)
c 
c       with the matrices b, v having dimensionalities
c       b(n,ncols), v(m,ncols), respectively. the reason for the
c       existence of this subroutine is the hope that the
c       dimensionality ncols, determined by this subroutine,
c       is comsiderably lower than either n or m
c 
c                     input parameters:
c 
c  a - the matrix to be decomposed (note that a is NOT damaged
c       by this subroutine in any way)
c  n,m - the dimensionalities of a
c  eps - the accuracy to which the decomposition (1) will be computed.
c 
c                     output parameters:
c 
c  b - the matrix in (1). note that the matrix b has to be dimensioned
c        b(n,m), but on exit from this subroutine only the first
c        ncols columns of  b  are meaningful, the rest making no
c        sense whatsoever, so that effectively on exit b is dimensioned
c        b(n,ncols)
c  v - the second matrix in (1). it is dimensioned v(m,ncols) on
c        exit. on entry the parameter ncols is unknown,
c        so v should be dimensioned v(m,n) by the user. however,
c        on exit only the first ncols columns of v are meaningful
c  ncols - the rank of the matrix a to the precision eps. also the
c        second dimension of the matrices b, v on exit
c  rnorms - the normalizing factors in the gram-schmidt process.
c        only the first ncols of them are meaninful, but the array
c        has to be dimensioned by the user to be at least m+1
c        real *8 elements long.
c 
c        . . . copy the user-supplied matrix a into b
c 
        do 1400 i=1,m
        do 1200 j=1,n
        b(j,i)=a(j,i)
 1200 continue
 1400 continue
c 
c        apply the gram-schmidt proces (with pivoting) to b
c 
         call rleasgrm(b,n,m,rnorms,eps,ncols,ifpivot)
c 
c        project the original matrix on the obtained orthogonal
c        columns
c 
        do 2400 i=1,ncols
        do 2200 j=1,m
c 
        call rleascap(a(1,j),b(1,i),n,prod)
cccc        v(j,i)=prod
        v(j,i)=(prod)
 2200 continue
 2400 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine rleasgrm(b,n,m,rnorms,eps,ncols,ifpivot)
        implicit real *8 (a-h,o-z)
        save
        real *8 b(n,m),cd
        dimension rnorms(1)
c 
c       this subroutine applies a pivoted double gram-schmidt
c       procedure to the matrix b, resulting in a collection of
c       orthogonal vectors. the number of these vectors is
c       equal to the numerical rank of b, given the precision
c       of computations eps, which is supplied by the user.
c 
c                    input paramneters:
c 
c  b - the matrix to be gram-schmidt'ed. it is destroyed by this
c       subroutine
c  n,m - dimensionalities of the matrix b
c  eps - the machine accuracy
c 
c                     output parameters:
c 
c  b - the matrix of gram-schmidt vectors of the matrix a. note
c        that on exit from this subroutine only the first
c        ncols columns of  b  are meaningful, the rest making no
c        sense whatsoever, so that effectively on exit b is dimensioned
c        b(n,ncols)
c  ncols - the rank of the matrix a to the precision eps. also the
c        second dimension of the matrix b on exit
c  rnorms - the normalizing factors in the gram-schmidt process.
c        only the first ncols of them are meaninful, but the array
c        has to be dimensioned by the user to be at least m+1
c        real *8 elements long.
c 
c       . . . prepare the array of values of norms
c             of columns
c 
        done=1
        dtot=0
        do 1400 i=1,m
        d=0
        do 1200 j=1,n
        d=d+b(j,i)*(b(j,i))
 1200 continue
        rnorms(i)=dsqrt(d)
        dtot=dtot+d
 1400 continue
c 
c       . . . conduct gram-schmidt iterations
c 
        thresh=dtot*eps**2
        do 4000 i=1,m
c 
c       find the pivot
c 
         if(ifpivot .eq. 0) goto 2700
        ipivot=i
        rn=rnorms(i)
c 
        do 2200 j=i+1,m
        if(rnorms(j) .le. rn) goto 2200
        rn=rnorms(j)
        ipivot=j
 2200 continue
 2400 continue
c 
c       put the column number ipivot in the i-th place
c 
        do 2600 j=1,n
        cd=b(j,i)
        b(j,i)=b(j,ipivot)
        b(j,ipivot)=cd
 2600 continue
c 
        d=rnorms(ipivot)
        rnorms(ipivot)=rnorms(i)
        rnorms(i)=d
 2700 continue
c 
c       orthogonalize the i-th column to all preceeding ones
c 
        if(i .eq. 1) goto 2790
        do 2780 j=1,i-1
c 
        call rleascap(b(1,i),b(1,j),n,cd)
c 
        do 2770 l=1,n
        b(l,i)=b(l,i)-b(l,j)*cd
 2770 continue
 2780 continue
 2790 continue
c 
c       normalize the i-th column
c 
        call rleascap(b(1,i),b(1,i),n,cd)
c 
c       if the norm of the remaining part of the matrix
c       is sufficiently small - terminate the process
c 
        d=cd
        if( (ifpivot .ne. 0) .and. (d .lt. thresh) ) return
        ncols=i
c 
        d=done/dsqrt(d)
        do 2800 j=1,n
        b(j,i)=b(j,i)*d
 2800 continue
c 
c        orthogonalize everything else to it
c 
        do 3200 j=i+1,m
c 
        if(rnorms(j) .lt. thresh/100) goto 3200
c 
        call rleascap(b(1,i),b(1,j),n,cd)
        cd=(cd)
c 
         rrn=0
        do 3000 l=1,n
        b(l,j)=b(l,j)-b(l,i)*cd
        rrn=rrn+b(l,j)*(b(l,j))
 3000 continue
        rnorms(j)=dsqrt(rrn)
 3200 continue
 3400 continue
 4000 continue
  
 4200 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine rleascap(x,y,n,prod)
        implicit real *8 (a-h,o-z)
        save
        dimension x(1),y(1)
c 
        prod=0
        do 1200 i=1,n
        prod=prod+x(i)*(y(i))
 1200 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine rleasrec(a,n,m,b)
        implicit real *8 (a-h,o-z)
        save
        real *8 a(n,m),b(n,m)
c 
c       reverse the columns of a
c 
        do 1400 i=1,m
        do 1200 j=1,n
        b(j,i)=a(j,m-i+1)
 1200 continue
 1400 continue
        return
c 
c 
c 
c 
        entry rleasrer(a,n,m,b)
c 
c       reverse the rows of a
c 
        do 2400 i=1,m
        do 2200 j=1,n
        b(j,i)=a(n-j+1,i)
 2200 continue
 2400 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine nrleamult(a,x,c,n,m,k)
        implicit real *8 (a-h,o-z)
        save
        dimension a(n,m),x(m,k),c(n,k)
c 
        do 1600 i=1,n
        do 1400 j=1,k
c 
        cd=0
        do 1200 ll=1,m
c 
        cd=cd+a(i,ll)*x(ll,j)
 1200 continue
c 
        c(i,j)=cd
 1400 continue
 1600 continue
        return
        end
  
c 
c 
c 
c 
        subroutine rleamat0(c,k,l,m,n,eps,delta,x,
     1      sa,ua,va,sb,ub,vb,work,lw,ltot,ncolsa,ncolsb)
        implicit real *8 (a-h,o-z)
c 
        save
        real *8 x(l,m),c(k,n),work(1),
     1      ua(k,1),va(l,1),ub(m,1),vb(n,1),sa(1),sb(1)
c 
c        this subroutine solves in the least squares sense the
c        equation
c 
c                     A(k,l) * X(l,m) * B(m,n) = C(k,n),               (1)
c 
c       where A, B, C are user-specified matrices, and X is the
c       matrix to be found. Note that the dimensionalities of
c       the matrices A,B,C, X are as general as could be.
c       Actually, this subroutine uses as its input the SVDs of
c       the matrices A,B, provided by the calling subroutine rleamatr
c       (see), so that the equation (1) has the form
c 
c                    UA * SA * VA^* * X * UB * SB * VB^* = C       (2)
c 
c       . . . multiply the matrix c from the left by UA^* ;
c             note that the resulting matrix UAC is dimensioned
c             uac(ncolsa,n)
c 
        call rleamatl(ier,ua,k,ncolsa,c,k,n,work)
c 
c       . . . multiply the matrix UAC by VB from the right;
c             note that the resulting matrix uacvb is
c             dimensioned uacvb(ncolsa,ncolsb)
c 
        call rleamat(ier,work,ncolsa,n,vb,n,ncolsb,x)
c 
c        now, the equation (2) has assumed the form
c 
c                   SA * VA^* * X * UB * SB  = UACVB              (3)
c 
c       . . . multiply (3) by SA^{-1} from the left and
c             by SB^{-1} from the right in the appropriate least
c             squares sense
c 
         call rleasas(sa,ncolsa,sb,ncolsb,x,eps,delta)
c 
c        now, the equation (3) has assumed the form
c 
c                  VA^* * X * UB  = UACVB                     (4)
c 
c       . . . multiply (4) by VA from the left and
c             by UB^* from the right, obtaining the
c             solution X
c 
        call rleamat(ier,va,l,ncolsa,x,ncolsa,ncolsb,work)
c 
        call rleamar(ier,work,l,ncolsb,ub,m,ncolsb,x)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine rleamat(ier,a,k,l,b,m,n,c)
        implicit real *8 (a-h,o-z)
        save
        real *8 a(k,l),b(m,n),c(k,n),d
c 
c        this subroutine multiplies the matrix a by the matrix b,
c        getting the matrix c
c 
c        . . . check that the dimensionalities of the input
c              matrices agree
c 
        ier=4
        if(l .eq. m) ier=0
        if(ier .ne. 0) return
c 
c       multiply the matrix a by the matrix b getting c
c 
        do 2000 i=1,k
        do 1800 j=1,n
        d=0
        do 1600 kk=1,l
        d=d+a(i,kk)*b(kk,j)
 1600 continue
        c(i,j)=d
 1800 continue
 2000 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine rleamatl(ier,a,k,l,b,m,n,c)
        implicit real *8 (a-h,o-z)
        save
        real *8 a(k,l),b(m,n),c(l,n),d
c 
c        this subroutine multiplies the matrix a^* by the matrix b,
c        getting the matrix c
c 
c        check that the dimensionalities of the input matrices agree
c 
        ier=4
        if(k .eq. m) ier=0
        if(ier .ne. 0) return
c 
c       multiply the adjoint of the matrix a by the matrix b getting c
c 
        do 2000 i=1,l
        do 1800 j=1,n
        d=0
        do 1600 kk=1,k
cccc        d=d+a(kk,i)*b(kk,j)
        d=d+(a(kk,i))*b(kk,j)
 1600 continue
        c(i,j)=d
 1800 continue
 2000 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine rleamar(ier,a,k,l,b,m,n,c)
        implicit real *8 (a-h,o-z)
        save
        real *8 a(k,l),b(m,n),c(k,m),d
c 
c        this subroutine multiplies the matrix a by the matrix b^*,
c        getting the matrix c
c 
c        check that the dimensionalities of the input matrices agree
c 
        ier=4
        if(l .eq. n) ier=0
        if(ier .ne. 0) return
c 
c       multiply the matrix a by the the adjoint of the matrix b getting c
c 
        do 2000 i=1,k
        do 1800 j=1,m
        d=0
        do 1600 kk=1,l
cccc        d=d+a(i,kk)*b(j,kk)
        d=d+a(i,kk)*(b(j,kk))
 1600 continue
        c(i,j)=d
 1800 continue
 2000 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine rleasas(sa,ncolsa,sb,ncolsb,uacvb,eps,delta)
        implicit real *8 (a-h,o-z)
        save
        real *8 sa(1),sb(1),uacvb(ncolsa,ncolsb)
c 
c       this subroutine multiplies the user-supplied matrix UACVB
c       by SA^{-1} from the left and by SB^{-1} from the right in
c       the appropriate least squares sense
c 
        do 1800 i=1,ncolsb
        do 1600 j=1,ncolsa
c 
        d=abs(sa(j)*sb(i))
        if(d .gt. eps) goto 1400
        uacvb(j,i)=0
        goto 1600
c 
 1400 continue
c 
        uacvb(j,i)=uacvb(j,i)/(delta+d)
 1600 continue
 1800 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine rleamem(x,y,n)
        implicit real *8 (a-h,o-z)
        save
        real *8 x(1),y(1)
c 
        do 1200 i=1,n
        y(i)=x(i)
 1200 continue
        return
        end
