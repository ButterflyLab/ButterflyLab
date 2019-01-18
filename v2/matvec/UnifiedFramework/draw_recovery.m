
%
% Reference: H. Yang, A Unified Framework for Oscillatory Integral 
% Transform: When to use NUFFT or Butterfly factorization? preprint, 2018.

close all;


%Set up parameters
N = 32; Nx = N; Np = N;
NG = 8;
NC = 256;
tol = 1e-12;
tR = 20; mR = 20; dim = 1;


k = -N/2:(N/2-1);
x = (0:N-1)/N;
kbox = [-N/2,N/2];
xbox = [0,1];
fun = @(x,k)fun_fio_1D(x,k);

kk = k(:);
xx = x(:);
K = fun(xx,kk);

%% test 1: a single butterfly matrix
Kt = K.';

% test the recovery of a complementary low-rank matrix via exact matvec
display('recovery of a complementary low-rank matrix via exact matvec');
% define matvec
fun = @(x) K*x;
funt = @(x) Kt*x;

% low-rank approximation of the phase function of an FIO
%  [U,S,V] = MVBF_Lowrank_Phase(N,N,fn,fnt,tol,tR,mR,dim,pi/2);

thre = pi/2;
isAL = 0;
type = 0;
% test discontinuity along row

rs = 1:floor(Nx/(tR+1)):Nx;%BF_RandSample(Nx,tR);
disPos2 = [];
for cnt = 1:tR
    tmp = zeros(Nx,1); tmp(rs(cnt)) = 1;
    data = funt(tmp);
    if isAL
        data = BF_analytic(data,2);
    end
    if type == 0
        data = data./abs(data);
    else
        data = conj(data)./abs(data);
    end
    M2 = imag(log(data));
    M2 = BF_Phase_Correction_Vec_1(M2,1,2*pi,thre);
    der2nd = M2(1:end-2)+M2(3:end)-2*M2(2:end-1);
    dis = find(abs(der2nd)>thre)+1;
    dis = dis(:)';
    disPos2 = unique([disPos2,dis]);
end

% test discontinuity along column
cs = 1:floor(Np/(tR+1)):Np;%BF_RandSample(Np,tR);
disPos1 = [];
for cnt = 1:tR
    tmp = zeros(Np,1); tmp(cs(cnt)) = 1;
    data = fun(tmp);
    if isAL
        data = BF_analytic(data,2);
    end
    data = data./abs(data);
    M1 = imag(log(data));
    M1 = BF_Phase_Correction_Vec_1(M1,1,2*pi,thre);
    der2nd = M1(1:end-2)+M1(3:end)-2*M1(2:end-1);
    dis = find(abs(der2nd)>thre)+1;
    dis = dis(:)';
    disPos1 = unique([disPos1,dis]);
end

disPos1 = unique([1,disPos1]);
disPos2 = unique([1,disPos2]);

rs = 1:floor(Nx/(tR+1)):Nx;%BF_RandSample(Nx,tR);
rs = [disPos1,disPos1+1,disPos1+2,rs]; rs = unique(rs); rs = sort(rs);
disPosSub1 = BF_find_entry(rs,disPos1);
tmp = zeros(Np,numel(rs));
tmp(rs+((1:numel(rs))-1)*Np) = 1;
data = funt(tmp);
if isAL
    data = BF_analytic(data,2);
end
if type == 0
    data = data./abs(data);
else
    data = conj(data)./abs(data);
end
M2 = imag(log(data));

%columns
cs = 1:floor(Np/(tR+1)):Np;%BF_RandSample(Np,tR);
cs = [disPos2,disPos2+1,disPos2+2,cs]; cs = unique(cs); cs = sort(cs);
disPosSub2 = BF_find_entry(cs,disPos2);
tmp = zeros(Nx,numel(cs));
tmp(cs+((1:numel(cs))-1)*Nx) = 1;
data = fun(tmp);
if isAL
    data = BF_analytic(data,2);
end
data = data./abs(data);
M1 = imag(log(data));

%[M1,M2] = BF_Phase_Correction_Mat_1(M1,M2,cs,rs,disPos1,disPos2,disPosSub1,disPosSub2,2*pi,dim);
U=M1;V=M2;posu=cs;posv=rs;disPosu=disPos1;disPosv=disPos2;disPosSubu=disPosSub1;disPosSubv=disPosSub2;tau=2*pi;
[mu,nu] = size(U);
[mv,nv] = size(V);

% correct first row in each continuous sector of rows
for cntd = 1:1%numel(disPosu)
    posc = disPosSubu(cntd);
    posc = 16;
    pic = figure;plot(V(:,posc),'-sb','LineWidth',2);axis square;
    axis square;axis([0,32,-5,40]);
    set(gca, 'FontSize', 32);
    b=get(gca);
    set(b.XLabel, 'FontSize', 32);set(b.YLabel, 'FontSize', 32);set(b.ZLabel, 'FontSize', 32);set(b.Title, 'FontSize', 32);
    tit = sprintf('./results/pic/p0.eps');
    saveas(pic,tit,'epsc');
    U = BF_Phase_Correction_Vec_1_draw(V(:,posc),disPosv,tau);
    pic = figure;plot(V(:,posc),'-sb','LineWidth',2);hold on ;plot(U,'-sr','LineWidth',2);axis square;
    axis square;axis([0,32,-5,40]);
    set(gca, 'FontSize', 32);
    b=get(gca);
    set(b.XLabel, 'FontSize', 32);set(b.YLabel, 'FontSize', 32);set(b.ZLabel, 'FontSize', 32);set(b.Title, 'FontSize', 32);
    tit = sprintf('./results/pic/p33.eps');
    saveas(pic,tit,'epsc');
end