close all;
clear all;

%Set up parameters
N = 1024;

func_name = 'funH';
switch func_name
    case 'funF' % Fourier transform
        k = -N/2:(N/2-1);
        x = (0:N-1)/N;
        kbox = [-N/2,N/2];
        xbox = [0,1];
        fun = @(x,k)funFT(x,k);
    case 'funIF' % Fourier transform
        k = -N/2:(N/2-1);
        x = (0:N-1)/N;
        kbox = [-N/2,N/2];
        xbox = [0,1];
        fun = @(x,k)funIFT(x,k);
    case 'fun0' % Fourier integral operator
        k = -N/2:(N/2-1);
        x = (0:N-1)/N;
        kbox = [-N/2,N/2];
        xbox = [0,1];
        fun = @(x,k)fun_fio_1D(x,k);
    case 'funH' % Hankel transformation
        k = 1:N;
        x = 1:N;
        kbox = [1,N+1];
        xbox = [1,N+1];
        fun = @(x,k)fun_Hankel(N,x,k);
end
kk = k(:);
xx = x(:);

%% Begin test
f = randn(N,1) + sqrt(-1)*randn(N,1);
tol = 1e-13;
mR = 10;

tic;
Factor = RSBFslow(fun, xx, xbox, kk, kbox, mR, tol, 1);
FactorT = toc;

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
disp(['Max Rank          : ' num2str(mR)]);
disp(['Tolerance         : ' num2str(tol)]);
disp(['Relative Error_2  : ' num2str(relerr)]);
disp(['Direct Time       : ' num2str(Td) ' s']);
disp(['Running Time      : ' num2str(RunT/60) ' mins']);
disp(['Factorization Time: ' num2str(FactorT/60) ' mins']);
disp(['Applying Time     : ' num2str(ApplyT) ' s']);
disp(['------------------------------------------']);

