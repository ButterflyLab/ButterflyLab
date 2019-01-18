function run_bf_explicit(N, func_name, mR, tol, fid)

addpath('../src/');
data_path = './data/';

if(~exist(data_path, 'dir'))
    mkdir(data_path);
end

kbox = [1,N+1];
k = 1:N;
kk = k(:);

xbox = [1,N+1];
x = 1:N;
xx = x(:);

switch func_name
    case 'fun0'
        fun = @(x,k)fun0(N,x,k);
    case 'funF'
        fun = @(x,k)funF(N,x,k);
    case 'funH'
        fun = @(x,k)funH(N,x,k);
end

f = randn(N,1) + 1i*randn(N,1);

tic;
Factor = bf_explicit(fun, xx, xbox, kk, kbox, mR, tol, 1);
FactorT = toc;

tic;
yy = apply_bf(Factor,f);
ApplyT = toc;
RunT = FactorT + ApplyT;

NC = 256;
tic;
relerr = bf_explicit_check(N,fun,f,yy,NC);
Td = toc;
Td = Td*N/NC;

fprintf(fid,'------------------------------------------\n');
fprintf(fid,'N                 : %4d\n', N);
fprintf(fid,'Max Rank          : %4d\n', mR);
fprintf(fid,'Tolerance         : %.3e\n', tol);
fprintf(fid,'Relative Error_2  : %.3e\n', relerr);
fprintf(fid,'Direct Time       : %.3e s\n', Td);
fprintf(fid,'Running Time      : %.3e mins\n', RunT/60);
fprintf(fid,'Factorization Time: %.3e mins\n', FactorT/60);
fprintf(fid,'Applying Time     : %.3e s\n', ApplyT);
fprintf(fid,'------------------------------------------\n\n');

save([data_path 'Factor_' func_name '_' num2str(N) '_' num2str(mR) '.mat'],'Factor','-v7.3');

end
