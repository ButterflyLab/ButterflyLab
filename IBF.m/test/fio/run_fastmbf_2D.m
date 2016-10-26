function run_fastmbf_2D(N, func_name, NG, tol, fid)

addpath('../../src/');
addpath('../kernels/');

k = -N/2:N/2-1;
[k1,k2] = ndgrid(k);
kk = [k1(:) k2(:)];

x = (0:N-1)/N;
[x1,x2] = ndgrid(x);
xx = [x1(:) x2(:)];

switch func_name
    case 'funFT'
        fun = @(x,k)funFT(x,k);
    case 'funIFT'
        fun = @(x,k)funIFT(x,k);
    case 'fun0'
        fun = @(x,k)fun0_2D(x,k);
    case 'fun1'
        fun = @(x,k)fun1_2D(x,k);
    case 'fun2'
        fun = @(x,k)fun2_2D(x,k);
end

f = randn(N^2,1) + 1i*randn(N^2,1);

tic;
[Factors,Rcomp] = fastMBF(fun,xx,kk,NG,tol);
FactorT = toc;

tic;
yy = apply_mbf(Factors,f);
ApplyT = toc;
RunT = FactorT + ApplyT;

NC = 256;
tic;
relerr = fbf_check(N,fun,f,xx,kk,yy,NC);
Td = toc;
Td = Td*N*N/NC;

fprintf(fid, '------------------------------------------\n');
fprintf(fid, 'N                 : %4d\n', N);
fprintf(fid, 'Chebyshev pts     : %4d\n', NG);
fprintf(fid, 'Tolerance         : %.3e\n', tol);
fprintf(fid, 'Relative Error_2  : %.3e\n', relerr);
fprintf(fid, 'Compression Ratio : %.3e\n', Rcomp);
fprintf(fid, 'Direct Time       : %.3e s\n', Td);
fprintf(fid, 'Running Time      : %.3e mins\n', RunT/60);
fprintf(fid, 'Factorization Time: %.3e mins\n', FactorT/60);
fprintf(fid, 'Applying Time     : %.3e s\n', ApplyT);
fprintf(fid, '------------------------------------------\n\n');

end
