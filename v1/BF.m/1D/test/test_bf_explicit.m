close all;
clear all;
clc;

addpath('../src/');
data_path = './data/';

if(~exist(data_path, 'dir'))
    mkdir(data_path);
end

N = 1024;
tol=1e-9;
mR = 4;

kbox = [1,N+1];
k = 1:N;
kk = k(:);

xbox = [1,N+1];
x = 1:N;
xx = x(:);

func_name = 'funF';
switch func_name
    case 'fun0'
        fun = @(x,k)fun0(N,x,k);
    case 'funF' % Fourier transform
        fun = @(x,k)funF(N,x,k);
    case 'funH' % Hankel function
        fun = @(x,k)funH(N,x,k);
end

%% Begin test
if(1)
    f = randn(N,1) + 1i*randn(N,1);
    binstr = sprintf('f_%d.bin', N);
    fid = fopen(binstr,'w');
    string = {'CpxNumMat'};
    serialize(fid, f, string);
end
if(0)
    binstr = sprintf('f_%d.bin', N);
    fid = fopen(binstr,'r');
    string = {'CpxNumMat'};
    f = deserialize(fid, string);
end

tic;
Factor = bf_explicit(fun, xx, xbox, kk, kbox, mR, tol, 1);
FactorT = toc;

tic;
yy = apply_bf(Factor,f);
ApplyT = toc;
RunT = FactorT + ApplyT;

NC = 256;
tic;
relerr = bf_explicit_check(N,fun,f,xx,kk,yy,NC);
Td = toc;
Td = Td*N/NC;

disp(['------------------------------------------']);
disp(['N                 : ' num2str(N)]);
disp(['Max Rank          : ' num2str(mR)]);
disp(['Tolerance         : ' num2str(tol)]);
disp(['Relative Error_2  : ' num2str(relerr)]);
disp(['Direct Time       : ' num2str(Td) ' s']);
disp(['Running Time      : ' num2str(RunT/60) ' mins']);
disp(['Factorization Time: ' num2str(FactorT/60) ' mins']);
disp(['Applying Time     : ' num2str(ApplyT) ' s']);
disp(['------------------------------------------']);

save([data_path 'Factor_' func_name '_' num2str(N) '_' num2str(mR) '.mat'],'Factor','-v7.3');

