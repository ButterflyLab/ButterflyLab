
%
% Reference: H. Yang, A Unified Framework for Oscillatory Integral 
% Transform: When to use NUFFT or Butterfly factorization? preprint, 2018.

%clc
close all;
clear all;

% This code compares two different approaches for the evaluation of 2D FIOs
% based on NUFFT


fid = fopen(['./ButterflyLab/results/Table_HankelNUFFT.log'],'a+');

for N = 2.^10%(8:2:14)
    NC = 256;
    tol = 1e-5;
    mR = 50; tR = 5*mR;
    
    k = (1:N)';
    x = (1:N)';
    kk = k(:);
    xx = x(:);
    fun = @(x,k)fun_Hankel(N,x,k);
    
    
    fn = @(y) fun(xx,kk(y));
    fnt = @(y) fun(xx(y),kk).';
    
    dim = 1;
    f = randn(N,1) + 1i*randn(N,1);
    rk = 2;
    numTest = 1;
    
    %% new NUFFT approach
    tic;
    % low-rank approximation and leading SVD of the phase function
    [U,S,V] = RSBF_Lowrank_Phase(N,N,fn,fnt,tol,mR,mR,dim,pi/2);
    V = V*S';
    recT = toc/60;
    
    % low rank approximation of the smooth part
    tic;
    funKsm = @(x,k) fun(xx(x,:),kk(k,:))./(exp(1i*(U(x,1:rk)*V(k,1:rk)')));
    [Us,Ss,Vs] = BF_lowrank(xx,kk,funKsm,tol,5*mR,mR);
    Vs = Ss*Vs';
    
    % use NUFFT to evaluate the matvec for exp(1i*U*V)
    if 0
        % greengard's NUFFT
        if rk == 1 % use 1D NUFFT
            nufftfun = @(cj) nufft1d3(N,V(:,1),cj,1,1e-12,N,U(:,1));
            oscFun = @(f) nufftfun(f);
        else if rk == 2 % use 2D NUFFT
                nufftfun = @(cj) nufft2d3(N,V(:,1),V(:,2),cj,1,1e-12,N,U(:,1),U(:,2));
                oscFun = @(f) nufftfun(f);
            else if rk == 3 % use 3D NUFFT
                    nufftfun = @(cj) nufft3d3(N,V(:,1),V(:,2),V(:,3),cj,1,1e-12,N,U(:,1),U(:,2),U(:,3));
                    oscFun = @(f) nufftfun(f);
                end
            end
        end
    else
        % use new NUFFT
        if rk == 1 % use 1D NUFFT
            nufftfun = nufftIII(V(:,1),U(:,1),1,N);
            oscFun = @(f) nufftfun(f);
        else if rk == 2 % use 2D NUFFT
                nufftfun = nufft2III(V(:,1:2),U(:,1:2),1,sqrt(N));
                oscFun = @(f) nufftfun(f)*N;
                 nufftfun2 = @(cj) nufft2d3(N,V(:,1),V(:,2),cj,1,1e-12,N,U(:,1),U(:,2));
                 oscFun2 = @(f) nufftfun2(f);
            else if rk == 3 % use 3D NUFFT
                    nufftfun = @(cj) nufft3d3(N,V(:,1),V(:,2),V(:,3),cj,1,1e-12,N,U(:,1),U(:,2),U(:,3));
                    oscFun = @(f) nufftfun(f);
                end
            end
        end
        
    end
    
    % apply the nonosc part
    rr = size(Us,2);
    maxRank = rr;
    appFun = @(f) sum(Us.*oscFun( (Vs.').*repmat(f,[1,rr]) ),2);
    FactorT = toc/60;
    
    appFun2 = @(f) sum(Us.*oscFun2( (Vs.').*repmat(f,[1,rr]) ),2);
    
    tic;
    for cnt = 1:numTest
        yy =  appFun(f);
    end
    ApplyT = toc/numTest;
%         yy2 =  appFun2(f);
    
    xs = BF_RandSample(N,NC);
    tic;
    yext = fun(xx(xs,:),kk)*f;
    Td = toc;
    Td = Td*N/NC;
    relerr = norm(yy(xs)-yext)/norm(yext);
%    relerr2 = norm(yy2(xs)-yext)/norm(yext)
%     yy(1:5)
%     yy2(1:5)

    disp(['------------------------------------------']);
    disp(['N                 : ' num2str(N)]);
    disp(['maxRank           : ' num2str(maxRank)]);
    disp(['Tolerance         : ' num2str(tol)]);
    disp(['Relative Error_2  : ' num2str(relerr)]);
    disp(['Direct Time       : ' num2str(Td) ' s']);
    disp(['Reconstruct Time  : ' num2str(recT) ' min']);
    disp(['Factorization Time: ' num2str(FactorT) ' min']);
    disp(['Applying Time     : ' num2str(ApplyT) ' s']);
    disp(['------------------------------------------']);
        
    
    fprintf(fid,'%7d &  %d &  %.2e & %.2e& %.2e & %.2e \\\\\n',...
        N,maxRank,recT,FactorT,ApplyT,relerr);
    
    
    
end
fprintf(fid,'\\ toprule \n');
fclose(fid);





