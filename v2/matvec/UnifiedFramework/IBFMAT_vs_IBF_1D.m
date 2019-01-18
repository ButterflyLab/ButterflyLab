% The code compares IBF-MAT and IBF.
%
% Copyright 2019 Haizhao Yang



close all;
clear all;

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
        fun = @(x,k)fun_fio_1D(x,k);
    case 'fun1'
        fun = @(x,k)fun_fio1_1D(x,k);
    case 'fun2'
        fun = @(x,k)fun_fio2_1D(x,k);
end

f = randn(N,1) + sqrt(-1)*randn(N,1);

tic;
[Factor,Rcomp] = IBF_Cheby(fun,xx,kk,NG,tol);
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
disp(['IBF']);
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

funk = @(t,s) fun(xx(t),kk(s));
tt = 1:N; tt = tt(:);
ss = 1:N; ss = ss(:);

tic;
[Factor,~] = IBF_uniform(funk,tt,ss,NG,tol);
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
disp(['IBF-MAT']);
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




