% This code tests the performance of the CUR butterfly factorizaiton for
% the 1D nonuniform Fourier transform.
%
% Copyright 2017 by Yingzhou Li and Haizhao Yang


% Set up parameters
i = 10;
N = 2^i;
tol = 1e-10;
NG = 8;  % number of Chebyshev pts

kk = (-N/2:N/2-1)';
%kk = rand(N,1)*(N-1)-N/2;

xx = ((0:N-1)/N)';
%xx = rand(N,1)*(N-1)/N;

fun = @(x,k)funFT(x,k);

f = randn(N,1) + sqrt(-1)*randn(N,1);

tic;
if 0
    [Factor,Rcomp] = RSBF(fun,xx,kk,NG,tol);
else
    funk = @(s,t) fun(xx(s),kk(t));
    ss = 1:N; ss = ss(:);
    tt = 1:N; tt = tt(:);
    [Factor,Rcomp] = IBF_uniform(funk,ss,tt,NG,tol);
end
FactorT = toc;

opCount = BF_op_count(Factor);

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
disp(['Chebyshev pts     : ' num2str(NG)]);
disp(['Tolerance         : ' num2str(tol)]);
disp(['Relative Error_2  : ' num2str(relerr)]);
disp(['Compression Ratio : ' num2str(Rcomp)]);
disp(['Prefactor         : ' num2str(opCount/N/log2(N))]);
disp(['Ratio wrt FFT     : ' num2str(opCount/N/log2(N)/34*9)]);
disp(['Direct Time       : ' num2str(Td) ' s']);
disp(['Running Time      : ' num2str(RunT/60) ' mins']);
disp(['Factorization Time: ' num2str(FactorT/60) ' mins']);
disp(['Applying Time     : ' num2str(ApplyT) ' s']);
disp(['------------------------------------------']);


data_path = './data/';
if(~exist(data_path, 'dir'))
    mkdir(data_path);
end
save([data_path 'Factor_nufft_' num2str(N) '_' num2str(NG) '_1D.mat'],'Factor','-v7.3');
