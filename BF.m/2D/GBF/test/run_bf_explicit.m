function run_bf_explicit(N, func_name, mR, tol, fid, saveflag)

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

switch func_name
    case 'funF'
        fun = @funF;
    case 'fun0'
        fun = @fun0;
    case 'fun1'
        fun = @fun1;
    case 'fun2'
        fun = @fun2;
end

f = randn(N,N) + 1i*randn(N,N);
f = reshape(f,N^2,1);

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
Td = Td*N^2/NC;

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

if(saveflag)
    save([data_path 'Factor_' func_name '_' num2str(N) '_' num2str(mR) '.mat'],'Factor','-v7.3');
end

end
