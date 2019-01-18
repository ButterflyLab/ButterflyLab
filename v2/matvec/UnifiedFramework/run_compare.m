clc
close all;
clear all;
%
% Reference: H. Yang, A Unified Framework for Oscillatory Integral 
% Transform: When to use NUFFT or Butterfly factorization? preprint, 2018.


%% Set up parameters
M = 512;
N = 1024;
tv = [1 4 16 64 256 1024];
[pp,pn,lamp,lamn] = waveEqn1D(M*8,M,N,tv);
NG = 12;
NC = 256;
tol = 1e-12;
mR = 50; tR = 5*mR; dim = 1;

f = randn(N,1);
numTest = 10;
tic;
for cnt = 1:numTest
    fft(f);
end
toc/numTest

k = -N/2:(N/2-1); kk = k(:);
x = (0:N-1)/N; xx = x(:);
for cntv = 1:numel(tv)
    
    % low-rank approximation of the phase function
    U = [pn{cntv},pp{cntv}]; V = [lamn,zeros(1,N/2);zeros(1,N/2),lamp];
    Ua = ones(N,1); Va = ones(1,N);
    
    
    %% IBF_uniform
    funk = @(t,s) exp(1i*U(t,:)*V(:,s));
    tt = 1:N; tt = tt(:);
    ss = 1:N; ss = ss(:);
    
    tic;
    [Factor,Rcomp] = IBF_uniform(funk,tt,ss,NG,tol);
    FactorT = toc;
    
    ampRk = size(Ua,2);
    f = randn(N,1) + sqrt(-1)*randn(N,1);
    tic;
    for cntt = 1:numTest
        yy = BF_apply(Factor,(Va.').*repmat(f,[1,ampRk]));
        yy = sum(Ua.*yy,2);
    end
    ApplyT = toc/numTest;
    RunT = FactorT + ApplyT;
    
    xs = BF_RandSample(N,NC);
    tic;
    yext = funk(xs,ss)*f;
    Td = toc;
    Td = Td*N/NC;
    relerr = norm(yy(xs)-yext)/norm(yext);
    
    disp(['------------------------------------------']);
    disp(['N                 : ' num2str(N)]);
    disp(['Interpolation pts     : ' num2str(NG)]);
    disp(['Tolerance         : ' num2str(tol)]);
    disp(['Relative Error_2  : ' num2str(relerr)]);
    disp(['Compress Rate     : ' num2str(Rcomp)]);
    disp(['Direct Time       : ' num2str(Td) ' s']);
    disp(['Running Time      : ' num2str(RunT) ' s']);
    disp(['Factorization Time: ' num2str(FactorT) ' s']);
    disp(['Applying Time     : ' num2str(ApplyT) ' s']);
    disp(['------------------------------------------']);
    
    %% NUFFT
    
    % compute the leading singular pairs of the phase function
    tic
    [Us,Ss,Vs] = BF_rsvd(U,V,1);
    Vs = Ss*Vs';
    
    rkThre = 50*ampRk; % determined by expected speed-up
    smp = BF_RandSample(N,rkThre*2);
    smp = sort(smp);
    amp = (Ua*Va(:,smp)).*exp(1i*U*V(:,smp))./exp(1i*Us*Vs(:,smp));
    ss = svd(amp,'econ');
    combRk = numel(find(ss>relerr*ss(1)));
    
    
    % decide whether we can use NUFFT
    if combRk/rkThre < 0.9
        isCan = 1;
    else
        isCan = 0;
    end
    % decide whether to use NUFFT
    if combRk < min(15*ampRk,rkThre) % determined by the comparison of NUFFT and BF
        isNUFFT = 1;
    else
        isNUFFT = 0;
    end
    timeDecide = toc;
    
    if 1%isCan
        if 0
            % visualize the smooth residual function
            smK = (Ua*Va).*exp(1i*U*V)./exp(1i*Us*Vs);
            figure;imagesc(real(smK));
        end
        tic;
        smKFun = @(x,k) (Ua(x,:)*Va(:,k)).*exp(1i*U(x,:)*V(:,k))./exp(1i*Us(x,:)*Vs(:,k));
        
        xidx = (1:N)'; kidx = xidx;
        [L,S,R] = BF_lowrank(xidx,kidx,smKFun,relerr,5*rkThre,rkThre);
        R = S*R';
        combRk = size(L,2);
        
        % use NUFFT to evaluate the matvec for exp(1i*Us*Vs)
        
        method = 1;
        switch method
            case 1 % greengard's NUFFT
                Uss = Us(:); Vss = Vs(:);
                nufftfun = @(cj) nufft1d3(N,Vss,cj,1,1e-12,N,Uss);
                ffun = @(f) sum(L.*nufftfun(((R.').*repmat(f,[1,combRk]))),2);
                timeFact = toc/60;
            case 2 % Alex's NUFFT
                sc = max(abs(Us));
                xxx = N*Us/sc;      %xxx = mod(xxx,N);
                kkk = Vs*sc/2/pi/N; kkk = kkk*N;%kkk = mod(kkk*N,N);
                nufftfun = nufftIII(-xxx(:),kkk(:),-1,N,mR,tol);
                ffun = @(f) sum(L.*nufftfun(((R.').*repmat(f,[1,ampRk]))),2);
                timeFact = toc/60;
            case 3
                sc = max(abs(Us));
                xxx = N*Us/sc;      %xxx = mod(xxx,N);
                kkk = Vs*sc/2/pi/N; kkk = kkk*N;%kkk = mod(kkk*N,N);
                ampFun = @(x,k) (Ua(x,:)*Va(:,k)).*exp(1i*U(x,:)*V(:,k))./exp(1i*Us(x,:)*Vs(:,k));
                [ffun,ampRk] = nufft1D3(ampFun,-xxx(:),kkk(:),-1,N,mR,relerr);
                timeFact = toc/60;
        end        % check new kernel
        if 0
            Knew = exp(2*pi*1i*xxx(:)*kkk(:)'/N);
            yy = Knew * f;
            norm(Knew-exp(1i*Us*Vs))/norm(exp(1i*Us*Vs))
        end
        
        tic;
        for cntt = 1:numTest
            y = ffun(f);
        end
        timeApply = toc/numTest;
        err = norm(y(xs)-yext)/norm(yext);
        
    else
        log10(ss')-log10(ss(1))
    end
    disp(['------------------------------------------']);
    disp(['N                 : ' num2str(N)]);
    disp(['isCan             : ' num2str(isCan)]);
    disp(['isNUFFT           : ' num2str(isNUFFT)]);
    disp(['combRk            : ' num2str(combRk)]);
    disp(['Relative Error_2  : ' num2str(err)]);
    disp(['Decision Time     : ' num2str(timeDecide) ' s']);
    disp(['Factorization Time: ' num2str(timeFact) ' s']);
    disp(['Applying Time     : ' num2str(timeApply) ' s']);
    disp(['------------------------------------------']);
end






