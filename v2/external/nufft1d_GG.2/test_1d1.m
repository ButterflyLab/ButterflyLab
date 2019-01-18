nj = 10000;
ms = 10000;

xj = (0:1/nj:(1-1/nj))*2*pi;
cj = exp(3*1i*xj);

eps=1e-12;


iflag = +1;

tic
fk = dirft1d1(nj,xj,cj,iflag,ms);
toc

tic
fk1 = nufft1d1(nj,xj,cj,iflag,eps,ms);
toc

abs_error=norm(fk-fk1,2)
rel_error=norm(fk-fk1,2)/norm(fk,2)


iflag = 1;

tic
fk = dirft1d1(nj,xj,cj,iflag,ms);
toc

tic
fk1 = nufft1d1(nj,xj,cj,iflag,eps,ms);
toc
fk1 = fk1([1,end:-1:2]);

tic;ff=fft(cj);toc
figure;plot(real(fftshift(ff(:)))/nj-real(fk1(:)));

abs_error=norm(fk-fk1,2)
rel_error=norm(fk-fk1,2)/norm(fk,2)
