function run_mbf_implicit(N, fun, fun_adj, mR, tol, fid)

addpath('../../GBF/src/');
addpath('../src/');
data_path = './data/';

if(~exist(data_path, 'dir'))
    mkdir(data_path);
end

kbox = [-N/2,N/2;-N/2,N/2];
k = -N/2:N/2-1;
[k1s,k2s] = ndgrid(k);
k1s = k1s(:);  k2s = k2s(:);
kk = [k1s k2s];

xbox = [0,1;0,1];
x = (0:N-1)'/N;
[x1s,x2s] = ndgrid(x);
x1s = x1s(:);  x2s = x2s(:);
xx = [x1s x2s];

f = randn(N,N) + 1i*randn(N,N);
f = reshape(f,N^2,1);

y = fun(f);

[Factor,FactorT] = mbf_implicit(fun, fun_adj, xx, xbox, kk, kbox, mR, tol, 0, 5);
if(FactorT<0)
    return;
end

tic;
yy = apply_mbf(Factor,f);
ApplyT = toc;
RunT = FactorT + ApplyT;

NC = 256;
relerr = mbf_implicit_check(N,y,yy,NC);

fprintf(fid,'------------------------------------------\n');
fprintf(fid,'N                 : %4d\n', N);
fprintf(fid,'Max Rank          : %4d\n', mR);
fprintf(fid,'Tolerance         : %.3e\n', tol);
fprintf(fid,'Relative Error_2  : %.3e\n', relerr);
fprintf(fid,'Running Time      : %.3e mins\n', RunT/60);
fprintf(fid,'Factorization Time: %.3e mins\n', FactorT/60);
fprintf(fid,'Applying Time     : %.3e s\n', ApplyT);
fprintf(fid,'------------------------------------------\n\n');

end
