% The IBF is credited to Yingzhou Li and Haizhao Yang, Interpolative 
% Butterfly Factorization, SIAM J. Sci. Comput., 39(2), A503-A531. 
%
% The IBF is a refined factorization version of the algorithm in E. J. 
% Candes, L. Demanet, and L. Ying. A fast butterfly algorithm for the 
% computation of Fourier integral operators. Multiscale Modeling and 
% Simulation, 7(4):1727?1750, 2009.
%
% Copyright 2017 by Yingzhou Li and Haizhao Yang



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
