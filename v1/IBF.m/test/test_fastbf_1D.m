close all;
clear all;

addpath('../src/');
addpath('./kernels/');

% Set up parameters
i = 10;
N = 2^i;
tol = 1e-6;
NG = 10;  % number of Chebyshev pts

k = -N/2:N/2-1;
kk = k(:);

x = (0:N-1)/N;
xx = x(:);

func_name = 'fun0';
switch func_name
    case 'funFT'
        fun = @(x,k)funFT(x,k);
    case 'funIFT'
        fun = @(x,k)funIFT(x,k);
    case 'fun0'
        fun = @(x,k)fun0_1D(x,k);
end

f = randn(N,1) + sqrt(-1)*randn(N,1);

tic;
[Factor,Rcomp] = fastBF(fun,xx,kk,NG,tol);
FactorT = toc;

tic;
yy = apply_fbf(Factor,f);
ApplyT = toc;
RunT = FactorT + ApplyT;

NC = 256;
tic;
relerr = fbf_check(N,fun,f,xx,kk,yy,NC);
Td = toc;
Td = Td*N/NC;

disp(['------------------------------------------']);
disp(['N                 : ' num2str(N)]);
disp(['Chebyshev pts     : ' num2str(NG)]);
disp(['Tolerance         : ' num2str(tol)]);
disp(['Relative Error_2  : ' num2str(relerr)]);
disp(['Compress Rate     : ' num2str(Rcomp)]);
disp(['Direct Time       : ' num2str(Td) ' s']);
disp(['Running Time      : ' num2str(RunT/60) ' mins']);
disp(['Factorization Time: ' num2str(FactorT/60) ' mins']);
disp(['Applying Time     : ' num2str(ApplyT) ' s']);
disp(['------------------------------------------']);


data_path = './data/';
if(~exist(data_path, 'dir'))
    mkdir(data_path);
end
save([data_path 'Factor_' func_name '_' num2str(N) '_' num2str(NG) '_1D.mat'],'Factor','-v7.3');
