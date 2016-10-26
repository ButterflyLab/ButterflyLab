close all;
clear all;
clc;

addpath('../src/');
data_path = './data/';

if(~exist(data_path, 'dir'))
    mkdir(data_path);
end

%% Set up parameters
N = 64;
tol=1e-4;
mR = 8;

% The index is different from 1D case.
% Here the range is from -N/2 to N/2-1 for each dimension of k
% and 0 to 1-1/N for each dimension of x.

kbox = [-N/2,N/2;-N/2,N/2];
k = -N/2:N/2-1;
[k1s,k2s] = ndgrid(k);
k1s = k1s(:);  k2s = k2s(:);
kk = [k1s k2s];

xbox = [0,1;0,1];
x = (0:N-1)'/N;
[x1s,x2s] = ndgrid(x);
x1s = x1s(:);  x2s = x2s(:);
xx = [x1s x2s];

func_name = 'funF';
switch func_name
    case 'funF'
        fun = @funF;
    case 'fun0'
        fun = @fun0;
    case 'fun1'
        fun = @fun1;
    case 'fun2'
        fun = @fun2;
end

%% Begin test
if(1)
    f = randn(N,N) + 1i*randn(N,N);
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
f = reshape(f,N^2,1);

tic;
Factor = bf_explicit(fun, xx, xbox, kk, kbox, mR, tol, 1);
FactorT = toc;

tic;
yy = apply_bf(Factor,f);
ApplyT = toc;
RunT = FactorT + ApplyT;

NC = 256;
tic;
relerr = bf_explicit_check(N,fun,f,yy,NC);
Td = toc;
Td = Td*N^2/NC;

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
