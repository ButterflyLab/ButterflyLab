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
sk = (0:n-1)';
tk = (0:n-1)';
uk = (0:n-1)';

eps=1e-12;


iflag = +1;

tic
fk = dirft3d3(nj,xj,yj,zj,cj,iflag,nk,sk,tk,uk);
toc

tic
fk1 = nufft3d3(nj,xj,yj,zj,cj,iflag,eps,nk,sk,tk,uk);
toc

abs_error=norm(fk-fk1,2)
rel_error=norm(fk-fk1,2)/norm(fk,2)


iflag = -1;

% tic
% fk = dirft3d3(nj,xj,yj,zj,cj,iflag,nk,sk,tk,uk);
% toc

tic
fk1 = nufft3d3(nj,xj,yj,zj,cj,iflag,eps,nk,sk,tk,uk);
toc

tmp = randn(20,20,20);
tic
ff = fftn(tmp);
toc

abs_error=norm(fk-fk1,2)
rel_error=norm(fk-fk1,2)/norm(fk,2)
