% The evaluation of Schlomilch expansions
%           f_k = sum_{n=1}^N c_n J_v(r_k w_n),  1<=k<=N
% where w_n = (n+v)pi, and 0<=r_1<...<r_N<=1
%
% The CURBF is credited to Eric Michielssen and Amir Boag, MULTILEVEL 
% EVALUATION OF ELECTROMAGNETIC FIELDS FOR THE RAPID SOLUTION OF SCAlTERlNG
% PROBLEMS, Microwave Opt. Technol. Lett., vol. 7, pp. 790-795, Dec. 1994.
%
% Copyright 2017 by Yingzhou Li and Haizhao Yang



function test_CUR_Schlomilch ( )

% Set up parameters
func_name = 'Schlomilch';
nu = 0;
i = 10;
N = 2^i;
tol = 1e-10;
NG = 10;  % number of Chebyshev pts
nv = 0;

w = ((1:N)+nv)*pi;
ww = w(:);

r = (0:N-1)/N;
rr = r(:);

f = randn(N,1) + sqrt(-1)*randn(N,1);

fun = @(r,w) fun_Bessel(nu,r,w);

tic;
[Factor,Rcomp] = CURBF(fun,rr,ww,NG,tol);
FactorT = toc;

tic;
yy = BF_apply(Factor,f);
ApplyT = toc;
RunT = FactorT + ApplyT;

NC = 256;
tic;
relerr = BF_check(N,fun,f,rr,ww,yy,NC);
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
return
end
