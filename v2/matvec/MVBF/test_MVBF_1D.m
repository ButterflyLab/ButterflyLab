close all;
clear all;

%Set up parameters
N = 1024;
NG = 8;
NC = 256;
tol = 1e-12;
tR = 20; mR = 20; dim = 1;

func_name = 'fun2';
switch func_name
    case 'funN' 
        c = 0.001;
        k = -N/2:(N/2-1);
        x = (0:N-1)/N;
        kbox = [-N/2,N/2];
        xbox = [0,1];
        fun = @(x,k)fun_fio3_1D(x,k,c);
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
if 1
    Kt = K.';
    
    % test the recovery of a complementary low-rank matrix via exact matvec
    display('recovery of a complementary low-rank matrix via exact matvec');
    % define matvec
    fn = @(x) K*x;
    fnt = @(x) Kt*x;
    
    % low-rank approximation of the phase function of an FIO
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
    
    % test the recovery of a complementary low-rank matrix via BF matvec
    display('recovery of a complementary low-rank matrix via BF matvec');
    
    if 1
        funK = @(x,k) K(x,k);
        [Factor,Rcomp] = IBF_uniform(funK,(1:N)',(1:N)',NG,tol);
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
    method = 1;
    switch method
        case 1 % IBF_uniform
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
        case 2 % CURBF
            V = S*V'; Va = Sa*Va';
            funk = @(t,s) exp(1i*U(t,:)*V(:,s));
            tt = 1:N; tt = tt(:);
            ss = 1:N; ss = ss(:);
            
            tic;
            [Factor,Rcomp] = CURBF(funk,tt,ss,NG,tol);
            FactorT = toc;
            
            rk = size(Ua,2);
            tic;
            yy = BF_apply(Factor,(Va.').*repmat(f,[1,rk]));
            yy = sum(Ua.*yy,2);
            ApplyT = toc;
            RunT = FactorT + ApplyT;
        case 3 % reference with an explicit kernel function
            tic;
            [Factor2,Rcomp] = IBF(fun,xx,kk,NG,tol);
            FactorT = toc;
            
            tic;
            yy = BF_apply(Factor2,f);
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
end

%% test 2: a composition of two FIOs
if func_name == 'fun2'
    display('Test2: a composition of two FIOs');
    % construct fast matvec
    % Factor = RSBFslow(fun, xx, xbox, kk, kbox, mR, tol, 0);
    if 1
        funK = @(x,k) K(x,k);
        [Factor0,Rcomp] = IBF_uniform(funK,(1:N)',(1:N)',NG,tol);
    else
        fun0 = @(t,s) K(t,s);
        Factor0 = CURBF(fun0,(1:N)',(1:N)',NG,tol);
    end
    
    fioComp = @(x)BF_apply(Factor0,fft(BF_apply(Factor0,x))/sqrt(N));
    fioCompt = @(y)BF_adj_apply(Factor0,ifft(BF_adj_apply(Factor0,y))*sqrt(N));
    K = K*fft(K)/sqrt(N);
    
    % low-rank approximation of the phase function
    [U,S,V] = MVBF_Lowrank_Phase(N,N,fioComp,fioCompt,tol,tR,mR,dim,pi/2,1);
    
    % low-rank approximation of the amplitude function
    [Ua,Sa,Va] = MVBF_Lowrank_Amp(N,N,fioComp,fioCompt,tol,tR,mR,dim);
    
    % check accuracy
    display('accuracy of the low-rank recovery:');
    display('overall accuracy:');
    norm((Ua*Sa*Va').*exp(1i*U*S*V')-K)
    display('phase accuracy:');
    norm(exp(1i*U*S*V')-K./abs(K))
    display('amplitude accuracy:');
    norm((Ua*Sa*Va')-abs(K))
    
    %% optimal BF via exact matvec and IBF or CURBF
    display('optimal BF via random sampling and IBF');
    f = randn(N,1) + sqrt(-1)*randn(N,1);
    method = 1;
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
            
        case 2 % IBF_uniform with a simple partition dealing with the discontinuity of FIOs, but the gap grows with N, otherwise error is large
            V = S*V'; Va = Sa*Va';
            gap = max(NG,32);
            par1 = 1:(N/2-gap);
            par2 = (N/2+1+gap):N;
            par3 = (N/2-gap+1):(N/2+gap);
            funk1 = @(t,s) exp(1i*U(t,:)*V(:,s));
            funk2 = @(t,s) exp(1i*U(t,:)*V(:,s+N/2+gap));
            tt = 1:N; tt = tt(:);
            ss = 1:(N/2-gap); ss = ss(:);
            
            tic;
            [Factor1,~] = IBF_uniform(funk1,tt,ss,NG,tol);
            [Factor2,~] = IBF_uniform(funk2,tt,ss,NG,tol);
            FactorT = toc;
            
            rk = size(Ua,2);
            tic;
            yy1 = BF_apply(Factor1,(Va(:,par1).').*repmat(f(par1),[1,rk]));
            yy2 = BF_apply(Factor2,(Va(:,par2).').*repmat(f(par2),[1,rk]));
            yy3 = ((Ua*Va(:,par3)).*exp(1i*U*V(:,par3)))*f(par3);
            yy1 = sum(Ua.*yy1,2);
            yy2 = sum(Ua.*yy2,2);
            yy = yy1+yy2+yy3;
            ApplyT = toc;
            RunT = FactorT + ApplyT;
            
        case 3 % IBF_uniform ignoring the discontinuity of the FIO, larger error when NG is large
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
            
        case 4 % CURBF
            V = S*V'; Va = Sa*Va';
            funk = @(t,s) exp(1i*U(t,:)*V(:,s));
            tt = 1:N; tt = tt(:);
            ss = 1:N; ss = ss(:);
            
            tic;
            [Factor,Rcomp] = CURBF(funk,tt,ss,NG,tol);
            FactorT = toc;
            
            rk = size(Ua,2);
            tic;
            yy = BF_apply(Factor,(Va.').*repmat(f,[1,rk]));
            yy = sum(Ua.*yy,2);
            ApplyT = toc;
            RunT = FactorT + ApplyT;
    end
    
    tic;
    ext = fioComp(f);
    relerr = norm(ext-yy)/norm(ext);
    Td = toc;
    Td = Td*N/NC;
    
    disp(['------------------------------------------']);
    disp(['N                 : ' num2str(N)]);
    disp(['Interpolation pts     : ' num2str(NG)]);
    disp(['Tolerance         : ' num2str(tol)]);
    disp(['Relative Error_2  : ' num2str(relerr)]);
    disp(['Direct Time       : ' num2str(Td) ' s']);
    disp(['Running Time      : ' num2str(RunT/60) ' mins']);
    disp(['Factorization Time: ' num2str(FactorT/60) ' mins']);
    disp(['Applying Time     : ' num2str(ApplyT) ' s']);
    disp(['------------------------------------------']);
end


