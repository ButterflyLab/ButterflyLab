close all;
clear all;

%Set up parameters
N = 1024;

func_name = 'funF';
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
    case 'fun2' % Fourier integral operator
        k = -N/2:(N/2-1);
        x = (0:N-1)/N;
        kbox = [-N/2,N/2];
        xbox = [0,1];
        fun = @(x,k)fun_fio2_1D(x,k);
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
NG = 10;

% construct fast matvec
Factor = RSBFslow(fun, xx, xbox, kk, kbox, NG, tol, 0);

fun = @(x)BF_apply(Factor,fft(BF_apply(Factor,x))/sqrt(N));
fun_adj = @(y)BF_adj_apply(Factor,ifft(BF_adj_apply(Factor,y))*sqrt(N));

f = randn(N,1) + sqrt(-1)*randn(N,1);

% exact result
tic;
y = fun(f);
toc;

% approximate result
tic;
Factor = MVBFslow(fun, fun_adj, xx, xbox, kk, kbox, NG, tol, 0);
FactorT = toc;

tic;
yy = BF_apply(Factor,f);
ApplyT = toc;
RunT = FactorT + ApplyT;

relerr = norm(y-yy)/norm(y);

disp(['------------------------------------------']);
disp(['N                 : ' num2str(N)]);
disp(['Max Rank          : ' num2str(NG)]);
disp(['Tolerance         : ' num2str(tol)]);
disp(['Relative Error_2  : ' num2str(relerr)]);
disp(['Running Time      : ' num2str(RunT/60) ' mins']);
disp(['Factorization Time: ' num2str(FactorT/60) ' mins']);
disp(['Applying Time     : ' num2str(ApplyT) ' s']);
disp(['------------------------------------------']);
