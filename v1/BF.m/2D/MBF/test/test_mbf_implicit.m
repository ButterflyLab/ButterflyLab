close all;
clear all;
%clc;

addpath('../../GBF/src/');
addpath('../../GBF/test/');
addpath('../src/');
data_path = './data/';

if(~exist(data_path, 'dir'))
    mkdir(data_path);
end

%% Set up parameters
N = 64;
tol=1e-4;
mR = 22;

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

FmR = 28;
func1_name = 'fun0';
func2_name = 'fun0';
Factor1 = load([data_path 'Factor_' func1_name '_' num2str(N) '_' num2str(FmR) '.mat'],'Factor');
Factor2 = Factor1;%load([data_path 'Factor_' func2_name '_' num2str(N) '_' num2str(FmR) '.mat'],'Factor');

fun = @(x)apply_mbf(Factor2.Factor,...
    reshape(fft2(reshape(apply_mbf(Factor1.Factor,x),N,N,[]))/N,N^2,[]));
fun_adj = @(y)apply_mbf_adj(Factor1.Factor,...
    reshape(ifft2(reshape(apply_mbf_adj(Factor2.Factor,y),N,N,[]))*N,N^2,[]));

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

y = fun(f);

[Factor,FactorT] = mbf_implicit(fun, fun_adj, xx, xbox, kk, kbox, mR, tol, 1, 0.02);
if(FactorT<0)
    return;
end

tic;
yy = apply_mbf(Factor,f);
ApplyT = toc;
RunT = FactorT + ApplyT;

NC = 256;
relerr = mbf_implicit_check(N,y,yy,NC);

disp(['------------------------------------------']);
disp(['N                 : ' num2str(N)]);
disp(['Max Rank          : ' num2str(mR)]);
disp(['Tolerance         : ' num2str(tol)]);
disp(['Relative Error_2  : ' num2str(relerr)]);
disp(['Running Time      : ' num2str(RunT/60) ' mins']);
disp(['Factorization Time: ' num2str(FactorT/60) ' mins']);
disp(['Applying Time     : ' num2str(ApplyT) ' s']);
disp(['------------------------------------------']);

%save([data_path 'Factor_' func1_name '_' func2_name '_' num2str(N) '_' num2str(mR) '.mat'],'Factor','-v7.3');
