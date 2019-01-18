close all;
clear all;
clc;

addpath('../src/');
data_path = './data/';

if(~exist(data_path, 'dir'))
    mkdir(data_path);
end

%% Set up parameters
N = 1024;
tol=1e-5;
mR = 6;

kbox = [1,N+1];
k = 1:N;
kk = k(:);

xbox = [1,N+1];
x = 1:N;
xx = x(:);

FmR = 4;
func1_name = 'fun0';
func2_name = 'fun0';
Factor1 = load([data_path 'Factor_' func1_name '_' num2str(N) '_' num2str(FmR) '.mat'],'Factor');
Factor2 = load([data_path 'Factor_' func2_name '_' num2str(N) '_' num2str(FmR) '.mat'],'Factor');

fun = @(x)apply_bf(Factor2.Factor,fft(apply_bf(Factor1.Factor,x))/sqrt(N));
fun_adj = @(y)apply_bf_adj(Factor1.Factor,ifft(apply_bf_adj(Factor2.Factor,y))*sqrt(N));

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
y = fun(f);
toc;

tic;
Factor = bf_implicit(fun, fun_adj, xx, xbox, kk, kbox, mR, tol, 1);
FactorT = toc;

tic;
yy = apply_bf(Factor,f);
ApplyT = toc;
RunT = FactorT + ApplyT;

NC = 256;
relerr = bf_implicit_check(N,y,yy,NC);

disp(['------------------------------------------']);
disp(['N                 : ' num2str(N)]);
disp(['Max Rank          : ' num2str(mR)]);
disp(['Tolerance         : ' num2str(tol)]);
disp(['Relative Error_2  : ' num2str(relerr)]);
disp(['Running Time      : ' num2str(RunT/60) ' mins']);
disp(['Factorization Time: ' num2str(FactorT/60) ' mins']);
disp(['Applying Time     : ' num2str(ApplyT) ' s']);
disp(['------------------------------------------']);

save([data_path 'Factor_' func1_name '_' func2_name '_' num2str(N) '_' num2str(mR) '.mat'],'Factor','-v7.3');
