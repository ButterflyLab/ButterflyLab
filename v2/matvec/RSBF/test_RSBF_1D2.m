close all;
clear all;

%Set up parameters
N = 1024;%512*64;
NG = 8;

func_name = 'fun0';
switch func_name
    case 'funF' % Fourier transform
        k = -N/2:(N/2-1);
        x = (0:N-1)/N;
        fun = @(x,k)funFT(x,k);
    case 'funIF' % Fourier transform
        k = -N/2:(N/2-1);
        x = (0:N-1)/N;
        fun = @(x,k)funIFT(x,k);
    case 'fun0' % Fourier integral operator
        k = -N/2:(N/2-1);
        x = (0:N-1)/N;
        fun = @(x,k)fun_fio_1D(x,k);
    case 'funH' % Hankel transformation
        k = 1:N;
        x = 1:N;
        fun = @(x,k)fun_Hankel(N,x,k);
end
kk = k(:);
xx = x(:);

K = fun(xx,kk);
K = K*fft(eye(N))*K;
Kt = K.';


%% test the recovery of a complementary low-rank matrix via random sampling
display('recovery of a complementary low-rank matrix via random sampling');
% define sampling operators
fn = @(x) K(:,x);
fnt = @(x) K(x,:).';
tol = 1e-13;
tR = 20; mR = 20; dim = 1;


%% optimal BF via random sampling and IBF or CURBF
display('optimal BF via random sampling and IBF_uniform or CURBF');
f = randn(N,1) + sqrt(-1)*randn(N,1);
method = 3;
switch method
    case 1 % IBF_uniform without amplitude
        V = S*V'; Va = Sa*Va';
        funk = @(t,s) exp(1i*U(t,:)*V(:,s));
        tt = 1:N; tt = tt(:);
        ss = 1:N; ss = ss(:);
        
        tic;
        [Factor,Rcomp] = IBF_uniform(funk,tt,ss,NG,tol);
        FactorT = toc;
        
        rk = size(Ua,2);
        tic;
        yy = BF_apply(Factor,(Va.').*repmat(f,[1,rk]));
        yy = sum(Ua.*yy,2);
        ApplyT = toc;
        RunT = FactorT + ApplyT;
    case 2 % IBF_uniform with amplitude
        funk = @(t,s) K(t,s);
        tt = 1:N; tt = tt(:);
        ss = 1:N; ss = ss(:);
        
        tic;
        [Factor,Rcomp] = IBF_uniform(funk,tt,ss,NG,tol);
        FactorT = toc;
        
        tic;
        yy = BF_apply(Factor,f);
        ApplyT = toc;
        RunT = FactorT + ApplyT;
    case 3 % CURBF
        % note that CURBF uses more than NG columns and rows for obtain better low-rank factorization
       
        funk = @(t,s) K(t,s);
        tt = 1:N; tt = tt(:);
        ss = 1:N; ss = ss(:);
        
        tic;
        [Factor,Rcomp] = CURBF(funk,tt,ss,NG,tol);
        FactorT = toc;
        
        tic;
        yy = BF_apply(Factor,f);
        ApplyT = toc;
        RunT = FactorT + ApplyT;
    case 4 % reference with an explicit kernel function
        tic;
        [Factor2,Rcomp] = IBF_Cheby(fun,xx,kk,NG,tol);
        FactorT = toc;
        
        tic;
        yy = BF_apply(Factor2,f);
        ApplyT = toc;
        RunT = FactorT + ApplyT;
end

NC = 256;
tic;
relerr = BF_check(N,fun,f,xx,kk,yy,NC);
Td = toc;
Td = Td*N/NC;

disp(['------------------------------------------']);
disp(['N                 : ' num2str(N)]);
disp(['Interpolation pts : ' num2str(NG)]);
disp(['Tolerance         : ' num2str(tol)]);
disp(['Relative Error_2  : ' num2str(relerr)]);
disp(['Compress Rate     : ' num2str(Rcomp)]);
disp(['Direct Time       : ' num2str(Td) ' s']);
disp(['Running Time      : ' num2str(RunT/60) ' mins']);
disp(['Factorization Time: ' num2str(FactorT/60) ' mins']);
disp(['Applying Time     : ' num2str(ApplyT) ' s']);
disp(['------------------------------------------']);

