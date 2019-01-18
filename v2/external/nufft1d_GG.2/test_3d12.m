n = 20;
nj = 20^3;
nk = 20^3;

xj = (0:1/n:(1-1/n))*2*pi;
yj = (0:1/n:(1-1/n))*2*pi;
zj = (0:1/n:(1-1/n))*2*pi;
[xj,yj,zj] = ndgrid(xj,yj,zj);
cj = exp(4*1i*(xj+yj+zj));
cc = cj;
xj = xj(:); yj = yj(:); zj = zj(:);
cj = cj(:);
ms = n;
mt = n;
mu = n;

eps=1e-12;


iflag = +1;

tic
fk = dirft3d1(nj,xj,yj,zj,cj,iflag,ms,mt,mu);
toc

tic
fk1 = nufft3d1(nj,xj,yj,zj,cj,iflag,eps,ms,mt,mu);
toc

abs_error=norm(fk-fk1,2)
rel_error=norm(fk-fk1,2)/norm(fk,2)


iflag = -1;

tic
fk = dirft3d1(nj,xj,yj,zj,cj,iflag,ms,mt,mu);
toc

tic
fk1 = nufft3d1(nj,xj,yj,zj,cj,iflag,eps,ms,mt,mu);
t1 = toc

abs_error=norm(fk-fk1,2)
rel_error=norm(fk-fk1,2)/norm(fk,2)

tic
ff = fftn(cc);
t2 = toc
figure;plot(real(fk1(:)));
figure;plot(fftshift(real(ff(:)))/nj);

rt = t1/t2