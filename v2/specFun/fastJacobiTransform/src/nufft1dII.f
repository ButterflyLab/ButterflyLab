	subroutine nufft1dII(nj,x,iflag,ns,rt,tol,U,V,xsub,r)
	implicit none

	integer :: ns,rt,iflag,nj,k,j,r,xsub(nj)
	real  :: tol
	real*8 pi,x(nj)
	parameter (pi=3.141592653589793238462643383279502884197d0)
	complex*16 fftconst,U(nj,r),V(ns,r)
	complex*16 M(nj,ns)

	fftconst = iflag*dcmplx(0,1)/ns*2*pi

	do k = 1,nj
	   do j = 1,ns
	      M(k,j) = exp(fftconst*(j-1)*(x(k)-floor(x(k)+0.5)))
	   enddo
	enddo

	!call lowrankfac(M,tol,rt,rt,U,V)
	xsub = mod(floor(x+0.5),ns)+1
	r = size(V,2)

	end subroutine

	subroutine nufft1dIIapp(nj,plan,c,U,V,xsub,ns,kflag,r,S)
	integer  r,i,j,k,nj,ns,kflag,num
	integer mm,xsub(nj)
	complex*16 M(r,ns),N(r,ns),S(nj),c(ns),c2(ns),U(r,nj),V(r,ns)
        complex*16,allocatable :: c1(:)
	double complex in1, out1
	dimension in1(ns), out1(nj)
	integer*8 :: plan

        c2=c
        if (kflag .lt. 0) then
           mm=floor(ns/2.0+0.6)
           allocate(c1(mm))
           c1=c2(1:mm)
           c2(1:mm)=c2(mm+1:ns)
           c2(mm+1:ns)=c1
        endif


	do k = 1,ns
	   M(:,k) = conjg(V(:,k))*c2(k)
	enddo
        !print *,'M(1,1:5)=',M(1,1:5)

	do i = 1,r
	   in1 = M(i,:)
	   call dfftw_execute_dft(plan, in1, out1)
	   N(i,:) = out1
	enddo
        !print *,'N(1,1:5)=',N(1,1:5)
        !print *,'xsub(1:5)',xsub(1:5)
	S = sum(U*N(:,xsub),1)

	end subroutine
