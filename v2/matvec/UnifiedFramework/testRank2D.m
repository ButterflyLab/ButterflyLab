
%
% Reference: H. Yang, A Unified Framework for Oscillatory Integral 
% Transform: When to use NUFFT or Butterfly factorization? preprint, 2018.

clear all;
close all;

i = 6;
N = 2^i;
tol = 1e-6;
NG = 8;  % number of Chebyshev pts

k = -N/2:N/2-1;
[k1,k2] = ndgrid(k);
kk = [k1(:) k2(:)];

x = (0:N-1)/N;
[x1,x2] = ndgrid(x);
xx = [x1(:) x2(:)];
fun = @(x,k)fun_fio_2D(x,k);
phFun = @(x,k)fun_phase_2D(x,k);
phFFT = @(x,k)fun_phase_FFT(x,k);

rd = sqrt(sum((kk').^2));
numW = 10;
st = 2*pi/numW; ed = 2*2*pi/numW;
pos = find(rd>0);
Pos_angle = zeros(N,N);
Pos_angle(pos) = acos(kk(pos,1)./rd(pos)');
        Pos_angle(:,1:end/2) = 2*pi - Pos_angle(:,1:end/2);
%figure;imagesc(Pos_angle);
idxMat = zeros(N,N);
pos = find(Pos_angle<ed & Pos_angle>st);
idxMat(pos) = 1;
idxMat(:,1:end/2) = 0;
%figure;imagesc(idxMat);
pos = find(idxMat(:)>0);

funK = @(x,k) fun(xx(x,:),kk(pos(k),:));
% do low-rank approximation
funPha = @(x,k) phFun(xx(x,:),kk(pos(k),:));
funPhaFT = @(x,k) phFFT(xx(x,:),kk(pos(k),:));
ix = (1:N^2)'; ip = (1:numel(pos))';
mR = 50;
[Uph,Sph,Vph] = BF_lowrank(ix,ip,funPha,1e-12,5*mR,mR);
Vph = Sph*Vph';
[U,S,V] = BF_rsvd(Uph,Vph,10);
ss = diag(S);
ss'
[Uph2,Sph2,Vph2] = BF_lowrank(ix,ip,funPhaFT,1e-12,5*mR,mR);
Vph2 = Sph2*Vph2';
[U2,S2,V2] = BF_rsvd(Uph2,Vph2,10);
ss2 = diag(S2);
ss2'

rk = 2;
funKsm = @(x,k) fun(xx(x,:),kk(pos(k),:))./(exp(1i*(U(x,1:rk)*S(1:rk,1:rk)*V(k,1:rk)')));
[Us,Ss,Vs] = BF_lowrank(ix,ip,funKsm,1e-12,5*mR,mR);
Vs = Ss*Vs';
[U2,S2,V2] = BF_rsvd(Us,Vs,50);
ss2 = diag(S2);
log10(ss2')


