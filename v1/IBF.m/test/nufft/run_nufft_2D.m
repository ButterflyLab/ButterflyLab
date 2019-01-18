function run_nufft_2D(N, NG, tol, fid)

addpath('../../src/');
addpath('../kernels/');

k = -N/2:N/2-1;
[k1,k2] = ndgrid(k);
kk = [k1(:) k2(:)];

if(~exist(sprintf('xx_%d_nufft_2D.bin', N), 'file'))
    fprintf('Generate non-uniform distribution of x from file\n');
    xx = rand(N^2,2)*(N-1)/N;
    binstr = sprintf('xx_%d_nufft_2D.bin', N);
    fidxx = fopen(binstr,'w');
    string = {'DblNumMat'};
    serialize(fidxx, xx, string);
    fclose(fidxx);
else
    fprintf('Read non-uniform distribution of x from file\n');
    binstr = sprintf('xx_%d_nufft_2D.bin', N);
    fidxx = fopen(binstr,'r');
    string = {'DblNumMat'};
    xx = deserialize(fidxx, string);
end

fun = @(x,k)funFT(x,k);

f = randn(N^2,1) + 1i*randn(N^2,1);

tic;
[Factor,Rcomp] = fastBF(fun,kk,xx,NG,tol);
FactorT = toc;

opCount = op_count(Factor);

tic;
yy = apply_fbf(Factor,f);
ApplyT = toc;
RunT = FactorT + ApplyT;

NC = 256;
tic;
relerr = fbf_check(N,fun,f,kk,xx,yy,NC);
Td = toc;
Td = Td*N*N/NC;

fprintf(fid, '------------------------------------------\n');
fprintf(fid, 'N                 : %4d\n', N);
fprintf(fid, 'Chebyshev pts     : %4d\n', NG);
fprintf(fid, 'Tolerance         : %.3e\n', tol);
fprintf(fid, 'Relative Error_2  : %.3e\n', relerr);
fprintf(fid, 'Compression Ratio : %.3e\n', Rcomp);
fprintf(fid, 'Prefactor         : %.3e\n', opCount/N^2/log2(N));
fprintf(fid, 'Ratio wrt FFT     : %.3e\n', opCount/N^2/log2(N)/34*9/2);
fprintf(fid, 'Direct Time       : %.3e s\n', Td);
fprintf(fid, 'Running Time      : %.3e mins\n', RunT/60);
fprintf(fid, 'Factorization Time: %.3e mins\n', FactorT/60);
fprintf(fid, 'Applying Time     : %.3e s\n', ApplyT);
fprintf(fid, '------------------------------------------\n\n');

end
