clc
close all;
clear all;

% This code compare two different FIO evaluation method: 1) BF; 2) NUFFT.
%
% Reference: H. Yang, A Unified Framework for Oscillatory Integral 
% Transform: When to use NUFFT or Butterfly factorization? preprint, 2018.

%% Set up FIO's of the wave equation
M = 512;
N = 1024;
tv = [1 16 256 4096];
[pp,pn,lamp,lamn] = waveEqn1D(M*8,M,N,tv);
NG = 12;
NC = 256;
tol = 1e-12;
mR = 50; tR = 5*mR;

f = randn(N,1);
numTest = 10;
tic;
for cnt = 1:numTest
    fft(f);
end
timeFFT = toc/numTest;

for cntv = 1:numel(tv)
    
    % low-rank approximation of the phase function
    U = [pn{cntv},pp{cntv}]; V = [lamn,zeros(1,N/2);zeros(1,N/2),lamp];
    
    
    %% IBF_uniform
    funk = @(t,s) exp(1i*U(t,:)*V(:,s));
    tt = 1:N; tt = tt(:);
    ss = 1:N; ss = ss(:);
    
    tic;
    [Factor,Rcomp] = IBF_uniform(funk,tt,ss,NG,tol);
    FactorT = toc;
    
    f = randn(N,1) + sqrt(-1)*randn(N,1);
    tic;
    for cntt = 1:numTest
        yy = BF_apply(Factor,f);
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
    
    %% NUFFT approach is based on the observation that the phase matrix is block-wise rank-1
    numBlk = 2;
    ffun = cell(1,numBlk);
    timeFact = 0;
    timeDecide = 0;
    for cntb = 1:numBlk
        % compute the leading singular pairs of the phase function
        vIdx = (cntb-1)*N/numBlk+(1:N/numBlk);
        tic
        [Us,Ss,Vs] = BF_rsvd(U,V(:,vIdx),1);
        Vs = Ss*Vs';
        
        rkThre = 50; % determined by expected speed-up
        smp = BF_RandSample(N/numBlk,rkThre*2);
        smp = sort(smp);
        amp = exp(1i*U*V(:,vIdx(smp)))./exp(1i*Us*Vs(:,smp));
        ss = svd(amp,'econ');
        combRk = numel(find(ss>relerr*ss(1)));
        
        
        % decide whether we can use NUFFT
        if combRk/rkThre < 0.9
            isCan = 1;
        else
            isCan = 0;
        end
        timeDecide = timeDecide + toc;
        
        tic;
        % use NUFFT to evaluate the matvec for exp(1i*Us*Vs)
        % greengard's NUFFT
        if size(Us,2) == 1 % use 1D NUFFT
            Uss = Us(:); Vss = Vs(:);
            nufftfun = @(cj) nufft1d3(N/numBlk,Vss,cj,1,1e-12,N,Uss);
            ffun{cntb} = @(f) nufftfun(f);
            timeFact = timeFact + toc;
        else if size(Us,2) == 2 % use 2D NUFFT
                nufftfun = @(cj) nufft2d3(N/numBlk,Vs(1,:)',Vs(2,:)',cj,1,1e-12,N,Us(:,1),Us(:,2));
                ffun{cntb} = @(f) nufftfun(f);
                timeFact = timeFact + toc;
            else if size(Us,2) == 3 % use 3D NUFFT
                nufftfun = @(cj) nufft3d3(N/numBlk,Vs(1,:)',Vs(2,:)',Vs(3,:)',cj,1,1e-12,N,Us(:,1),Us(:,2),Us(:,3));
                ffun{cntb} = @(f) nufftfun(f);
                timeFact = timeFact + toc;
                end
            end
        end
    end
    
    tic;
    for cntt = 1:numTest
        y = zeros(N,1);
        for cntb = 1:numBlk
            vIdx = (cntb-1)*N/numBlk+(1:N/numBlk);
            y = y + ffun{cntb}(f(vIdx));
        end
    end
    timeApply = toc/numTest;
    err = norm(y(xs)-yext)/norm(yext);
    
    disp(['------------------------------------------']);
    disp(['N                 : ' num2str(N)]);
    disp(['isCan             : ' num2str(isCan)]);
    disp(['combRk            : ' num2str(combRk)]);
    disp(['Relative Error_2  : ' num2str(err)]);
    disp(['FFT Time          : ' num2str(timeFFT) ' s']);
    disp(['Decision Time     : ' num2str(timeDecide) ' s']);
    disp(['Factorization Time: ' num2str(timeFact) ' s']);
    disp(['Applying Time     : ' num2str(timeApply) ' s']);
    disp(['------------------------------------------']);
end






