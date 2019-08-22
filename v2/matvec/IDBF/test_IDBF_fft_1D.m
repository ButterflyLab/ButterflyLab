% This code tests the performance of the IDBF for
% the 1D Fourier transform.
%
% Copyright 2019 Haizhao Yang


% Set up parameters
i = 10;
N = 2^i; % problem size
tol = 1e-10; % target approximation accuracy
rk = 16;  % rank parameter for low-rank matrices
n0 = 8; % number of grid points per leaf
tt = 5; % over sampling parameter in rand ID
rand_or_cheb = 'rand'; 
% 'rand' randomized sampling low-rank approximation
% 'cheb' mock-chebyshev low-rank approximation

% discretization of the FFT
kk = (-N/2:N/2-1)';
%kk = rand(N,1)*(N-1)-N/2;

xx = ((0:N-1)/N)';
%xx = rand(N,1)*(N-1)/N;

fun = @(x,k)funFT(x,k);

f = randn(N,1) + sqrt(-1)*randn(N,1);

% perform IDBF
tic;
Factor = BF_IDBF(fun,xx,kk,n0,rk,tol,rand_or_cheb,tt);
FactorT = toc;


% test the accuracy
tic;
yy = BF_apply(Factor,f);
ApplyT = toc;
RunT = FactorT + ApplyT;

NC = 256;
tic;
relerr = BF_check(N,fun,f,xx,kk,yy,NC);
Td = toc;
Td = Td*N/NC;

disp(['------------------------------------------']);
disp(['N                 : ' num2str(N)]);
disp(['Chebyshev pts     : ' num2str(rk)]);
disp(['Tolerance         : ' num2str(tol)]);
disp(['Relative Error_2  : ' num2str(relerr)]);
disp(['Direct Time       : ' num2str(Td) ' s']);
disp(['Running Time      : ' num2str(RunT/60) ' mins']);
disp(['Factorization Time: ' num2str(FactorT/60) ' mins']);
disp(['Applying Time     : ' num2str(ApplyT) ' s']);
disp(['------------------------------------------']);

