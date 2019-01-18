function run_MVBFslow(N, fun, fun_adj, NG, tol, fid)

%
% Reference: H. Yang, A Unified Framework for Oscillatory Integral 
% Transform: When to use NUFFT or Butterfly factorization? preprint, 2018.

k = -N/2:(N/2-1);
x = (0:N-1)/N;
kbox = [-N/2,N/2];
xbox = [0,1];
kk = k(:);
xx = x(:);

f = randn(N,1) + 1i*randn(N,1);

y = fun(f);

tic;
Factor = MVBFslow(fun, fun_adj, xx, xbox, kk, kbox, NG, tol, 0);
FactorT = toc;

tic;
yy = BF_apply(Factor,f);
ApplyT = toc;
RunT = FactorT + ApplyT;

relerr = norm(y-yy)/norm(y);

fprintf(fid,'------------------------------------------\n');
fprintf(fid,'N                 : %4d\n', N);
fprintf(fid,'Max Rank          : %4d\n', NG);
fprintf(fid,'Tolerance         : %.3e\n', tol);
fprintf(fid,'Relative Error_2  : %.3e\n', relerr);
fprintf(fid,'Running Time      : %.3e mins\n', RunT/60);
fprintf(fid,'Factorization Time: %.3e mins\n', FactorT/60);
fprintf(fid,'Applying Time     : %.3e s\n', ApplyT);
fprintf(fid,'------------------------------------------\n\n');

end
