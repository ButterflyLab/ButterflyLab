% This code tests the performance of the CUR butterfly factorizaiton for
% the 2D nonuniform Fourier transform.
%
% The CURBF is credited to Eric Michielssen and Amir Boag, MULTILEVEL 
% EVALUATION OF ELECTROMAGNETIC FIELDS FOR THE RAPID SOLUTION OF SCAlTERlNG
% PROBLEMS, Microwave Opt. Technol. Lett., vol. 7, pp. 790-795, Dec. 1994.
%
% Copyright 2017 by Yingzhou Li and Haizhao Yang


% Set up parameters
i = 6;
N = 2^i;
tol = 1e-8;
NG = 32;  % number of Chebyshev pts

k = -N/2:N/2-1;
[k1,k2] = ndgrid(k);
kk = [k1(:) k2(:)];

xx = rand(N^2,2)*(N-1)/N;

fun = @(x,k)funFT(x,k);

f = randn(N^2,1) + sqrt(-1)*randn(N^2,1);

tic;
[Factor,Rcomp] = CURBF(fun,xx,kk,NG,tol);
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
Td = Td*N*N/NC;

disp(['------------------------------------------']);
disp(['N                 : ' num2str(N)]);
disp(['Chebyshev pts     : ' num2str(NG)]);
disp(['Tolerance         : ' num2str(tol)]);
disp(['Relative Error_2  : ' num2str(relerr)]);
disp(['Compression Ratio : ' num2str(Rcomp)]);
disp(['Prefactor         : ' num2str(opCount/N^2/log2(N))]);
disp(['Ratio wrt FFT     : ' num2str(opCount/N^2/log2(N)/34*9/2)]);
disp(['Direct Time       : ' num2str(Td) ' s']);
disp(['Running Time      : ' num2str(RunT/60) ' mins']);
disp(['Factorization Time: ' num2str(FactorT/60) ' mins']);
disp(['Applying Time     : ' num2str(ApplyT) ' s']);
disp(['------------------------------------------']);


data_path = './data/';
if(~exist(data_path, 'dir'))
    mkdir(data_path);
end
save([data_path 'Factor_nufft_' num2str(N) '_' num2str(NG) '_2D.mat'],'Factor','-v7.3');
