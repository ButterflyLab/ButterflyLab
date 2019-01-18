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
i = 4;
N = 2^i;
tol = 1e-6;
NG = 8;  % number of Chebyshev pts

k = -N/2:N/2-1;
[k1,k2,k3] = ndgrid(k);
kk = [k1(:) k2(:) k3(:)];

xx = rand(N^3,3)*(N-1)/N;

fun = @(x,k)funFT(x,k);

f = randn(N^3,1) + 1i*randn(N^3,1);

tic;
[Factor,Rcomp] = IBF_Cheby(fun,kk,xx,NG,tol);
FactorT = toc;

opCount = BF_op_count(Factor);

tic;
yy = BF_apply(Factor,f);
ApplyT = toc;
RunT = FactorT + ApplyT;

NC = 256;
tic;
relerr = BF_check(N,fun,f,kk,xx,yy,NC);
Td = toc;
Td = Td*N*N*N/NC;

disp(['------------------------------------------']);
disp(['N                 : ' num2str(N)]);
disp(['Chebyshev pts     : ' num2str(NG)]);
disp(['Tolerance         : ' num2str(tol)]);
disp(['Relative Error_2  : ' num2str(relerr)]);
disp(['Compression Ratio : ' num2str(Rcomp)]);
disp(['Prefactor         : ' num2str(opCount/N^3/log2(N))]);
disp(['Ratio wrt FFT     : ' num2str(opCount/N^3/log2(N)/34*9/3)]);
disp(['Direct Time       : ' num2str(Td) ' s']);
disp(['Running Time      : ' num2str(RunT/60) ' mins']);
disp(['Factorization Time: ' num2str(FactorT/60) ' mins']);
disp(['Applying Time     : ' num2str(ApplyT) ' s']);
disp(['------------------------------------------']);


data_path = './data/';
if(~exist(data_path, 'dir'))
    mkdir(data_path);
end
save([data_path 'Factor_nufft_' num2str(N) '_' num2str(NG) '_3D.mat'],'Factor','-v7.3');
