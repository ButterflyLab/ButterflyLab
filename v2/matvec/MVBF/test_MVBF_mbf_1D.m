close all;
clear all;

%Set up parameters
N = 1024;
NG = 12;
NC = 256;
tol = 1e-12;
tR = 20; mR = 20; dim = 1;

func_name = 'fun2';
switch func_name
    case 'funF' % Fourier transform
        k = -N/2:(N/2-1);
        x = (0:N-1)/N;
        kbox = [-N/2,N/2];
        xbox = [0,1];
        fun = @(x,k)funFT(x,k);
    case 'funIF' % Fourier transform
        k = -N/2:(N/2-1);
        x = (0:N-1)/N;
        kbox = [-N/2,N/2];
        xbox = [0,1];
        fun = @(x,k)funIFT(x,k);
    case 'fun2' % Fourier integral operator
        k = -N/2:(N/2-1);
        x = (0:N-1)/N;
        kbox = [-N/2,N/2];
        xbox = [0,1];
        fun = @(x,k)fun_fio2_1D(x,k);
    case 'funH' % Hankel transformation
        k = 1:N;
        x = 1:N;
        kbox = [1,N+1];
        xbox = [1,N+1];
        fun = @(x,k)fun_Hankel(N,x,k);
end
kk = k(:);
xx = x(:);
K = fun(xx,kk);

%% test 1: a single butterfly matrix
% test the recovery of a complementary low-rank matrix via BF matvec
display('recovery of a complementary low-rank matrix via BF matvec');

if 1
    [Factor,Rcomp] = IBF_Cheby(fun,xx,kk,NG,tol);
else
    [Factor,Rcomp] = CURBF(fun,xx,kk,NG,tol);% TODO: find a better way for CUR
end

f = randn(N,1) + sqrt(-1)*randn(N,1);
yy = BF_apply(Factor,f);
relerr = BF_check(N,fun,f,xx,kk,yy,NC)

fn = @(f) BF_apply(Factor,f);
fnt = @(f) BF_trans_apply(Factor,f);

% low-rank approximation of the phase function of a matvec
[U,S,V] = MVBF_Lowrank_Phase(N,N,fn,fnt,tol,tR,mR,dim,pi/2);

% low-rank approximation of the amplitude function
[Ua,Sa,Va] = MVBF_Lowrank_Amp(N,N,fn,fnt,tol,tR,mR,dim);

% check accuracy
display('accuracy of the low-rank recovery:');
display('overall accuracy:');
norm((Ua*Sa*Va').*exp(1i*U*S*V')-K)
display('phase accuracy:');
norm(exp(1i*U*S*V')-K./abs(K))
display('amplitude accuracy:');
norm((Ua*Sa*Va')-abs(K))

%% optimal BF via exact matvec and IBF or CURBF
display('Test1: a single FIO or oscillatory matrix');
f = randn(N,1) + sqrt(-1)*randn(N,1);
method = 2;
switch method
    case 1 % IBF_multiscale_uniform
        V = S*V'; Va = Sa*Va';
        funk = @(t,s) exp(1i*U(t,:)*V(:,s));
        tt = 1:N; tt = tt(:);
        ss = 1:N; ss = ss(:);
        
        tic;
        [Factor,Rcomp] = IBF_multiscale_uniform(funk,tt,ss,NG,tol);
        FactorT = toc;
        
        rk = size(Ua,2);
        tic;
        yy = BF_multiscale_apply(Factor,(Va.').*repmat(f,[1,rk]));
        yy = sum(Ua.*yy,2);
        ApplyT = toc;
        RunT = FactorT + ApplyT;
    case 2 % IBF_multiscale_Cheby
        tic;
        [Factor,Rcomp] = IBF_multiscale_Cheby(fun,xx,kk,NG,tol);
        FactorT = toc;
        
        tic;
        yy = BF_multiscale_apply(Factor,f);
        ApplyT = toc;
        RunT = FactorT + ApplyT;
end

tic;
relerr = BF_check(N,fun,f,xx,kk,yy,NC);
Td = toc;
Td = Td*N/NC;

disp(['------------------------------------------']);
disp(['N                 : ' num2str(N)]);
disp(['Interpolation pts     : ' num2str(NG)]);
disp(['Tolerance         : ' num2str(tol)]);
disp(['Relative Error_2  : ' num2str(relerr)]);
disp(['Compress Rate     : ' num2str(Rcomp)]);
disp(['Direct Time       : ' num2str(Td) ' s']);
disp(['Running Time      : ' num2str(RunT/60) ' mins']);
disp(['Factorization Time: ' num2str(FactorT/60) ' mins']);
disp(['Applying Time     : ' num2str(ApplyT) ' s']);
disp(['------------------------------------------']);