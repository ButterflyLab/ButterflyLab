
%
% Reference: H. Yang, A Unified Framework for Oscillatory Integral 
% Transform: When to use NUFFT or Butterfly factorization? preprint, 2018.

clc
close all;
clear all;

% This code compares two different approaches for the evaluation of 2D FIOs
% based on NUFFT


fid = fopen(['./ButterflyLab/results/Table_FIONUFFT_2v2.log'],'a+');
fid2 = fopen(['./ButterflyLab/results/Table_FIONUFFT_2v3.log'],'a+');

for N = 2.^(4:10)
    %% Set up FIO's
    NC = 256;
    tol = 1e-12;
    mR = 100; tR = 5*mR;
    
    k = -N/2:N/2-1;
    [k1,k2] = ndgrid(k);
    kk = [k1(:) k2(:)];
    rd = sqrt(sum((kk').^2));
    
    pos = find(rd>0);
    Pos_angle = zeros(N,N);
    Pos_angle(pos) = acos(kk(pos,1)./rd(pos)');
    Pos_angle(:,1:end/2) = 2*pi - Pos_angle(:,1:end/2);
    
    x = (0:N-1)/N;
    [x1,x2] = ndgrid(x);
    xx = [x1(:) x2(:)];
    
    fun = @(x,k)fun_fio_2D(x,k);
    funp = @(x,k)fun_phase_2D(x,k);
    
    f = randn(N^2,1);
    numTest = 1;
    
    
    %% new NUFFT approach
    if 1
        rk = 2;
        numW = 8;
        oscFun = cell(1,numW);
        appFun = cell(1,numW);
        maxRank = 0;
        tic;
        for cntw = 1:numW
            dw = 2*pi/numW; st = (cntw-1)*dw; ed = st + dw;
            idxMat = zeros(N,N);
            pos = find(Pos_angle<ed & Pos_angle>=st);
            idxMat(pos) = 1;
            pos = find(idxMat(:)>0);
            wlen = numel(pos);
            
            % low rank approximation and leading SVD of the phase
            if 1
                funPha = @(x,k) funp(xx(x,:),kk(pos(k),:));
                ix = (1:N^2)'; ip = (1:numel(pos))';
                [Uph,Sph,Vph] = BF_lowrank(ix,ip,funPha,tol,5*mR,mR);
                if 0
                    Vph = Sph*Vph';
                    [U,S,V] = BF_rsvd(Uph,Vph,10);
                    U = U*S;
                else
                    U = Uph*sqrt(Sph); V = Vph*sqrt(Sph);
                end
            else
                agl = (st+ed)/2;
                kkc = [cos(agl),sin(agl)];
                Us = funp(xx,kkc);
                Vs = sqrt(sum(kk(pos,:).^2,2));
                
                funPha = @(x,k) funp(xx(x,:),kk(pos(k),:))-Us(x)*Vs(k)';
                ix = (1:N^2)'; ip = (1:numel(pos))';
                [Uph,Sph,Vph] = BF_lowrank(ix,ip,funPha,tol,5*mR,mR);
                U = Uph*Sph; V = Vph;
                U = [Us,U]; V = [Vs,V];
            end
            
            % low rank approximation of the smooth part
            funKsm = @(x,k) fun(xx(x,:),kk(pos(k),:))./(exp(1i*(U(x,1:rk)*V(k,1:rk)')));
            [Us,Ss,Vs] = BF_lowrank(ix,ip,funKsm,tol,5*mR,mR);
            Vs = Ss*Vs';
            
            % use NUFFT to evaluate the matvec for exp(1i*Us*Vs)
            % greengard's NUFFT
            if rk == 1 % use 1D NUFFT
                nufftfun = @(cj) nufft1d3(wlen,V(:,1),cj,1,1e-12,N^2,U(:,1));
                oscFun{cntw} = @(f) nufftfun(f);
            else if rk == 2 % use 2D NUFFT
                    nufftfun = @(cj) nufft2d3(wlen,V(:,1),V(:,2),cj,1,1e-12,N^2,U(:,1),U(:,2));
                    oscFun{cntw} = @(f) nufftfun(f);
                else if rk == 3 % use 3D NUFFT
                        nufftfun = @(cj) nufft3d3(wlen,V(:,1),V(:,2),V(:,3),cj,1,1e-12,N^2,U(:,1),U(:,2),U(:,3));
                        oscFun{cntw} = @(f) nufftfun(f);
                    end
                end
            end
            
            % apply the nonosc part
            rr = size(Us,2);
            maxRank = max(maxRank,rr);
            appFun{cntw} = @(f) sum(Us.*oscFun{cntw}( (Vs.').*repmat(f(pos),[1,rr]) ),2);
        end
        FactorT = toc;
        
        tic;
        for cnt = 1:numTest
            yy = zeros(N^2,1);
            for cntw = 1:numW
                yy = yy + appFun{cntw}(f);
            end
        end
        ApplyT = toc/numTest;
        
        xs = BF_RandSample(N^2,NC);
        tic;
        yext = fun(xx(xs,:),kk)*f;
        Td = toc;
        Td = Td*N/NC;
        relerr = norm(yy(xs)-yext)/norm(yext);
        
        disp(['------------------------------------------']);
        disp(['N                 : ' num2str(N)]);
        disp(['maxRank           : ' num2str(maxRank)]);
        disp(['Tolerance         : ' num2str(tol)]);
        disp(['Relative Error_2  : ' num2str(relerr)]);
        disp(['Direct Time       : ' num2str(Td) ' s']);
        disp(['Factorization Time: ' num2str(FactorT) ' s']);
        disp(['Applying Time     : ' num2str(ApplyT) ' s']);
        disp(['------------------------------------------']);
    end
    
    if 1
        rk = 3;
        numW = 8;
        oscFun = cell(1,numW);
        appFun = cell(1,numW);
        maxRank3 = 0;
        tic;
        for cntw = 1:numW
            dw = 2*pi/numW; st = (cntw-1)*dw; ed = st + dw;
            idxMat = zeros(N,N);
            pos = find(Pos_angle<ed & Pos_angle>=st);
            idxMat(pos) = 1;
            pos = find(idxMat(:)>0);
            wlen = numel(pos);
            
            % low rank approximation and leading SVD of the phase
            if 1
                funPha = @(x,k) funp(xx(x,:),kk(pos(k),:));
                ix = (1:N^2)'; ip = (1:numel(pos))';
                [Uph,Sph,Vph] = BF_lowrank(ix,ip,funPha,tol,5*mR,mR);
                if 0
                    Vph = Sph*Vph';
                    [U,S,V] = BF_rsvd(Uph,Vph,10);
                    U = U*S;
                else
                    U = Uph*sqrt(Sph); V = Vph*sqrt(Sph);
                end
            else
                agl = (st+ed)/2;
                kkc = [cos(agl),sin(agl)];
                Us = funp(xx,kkc);
                Vs = sqrt(sum(kk(pos,:).^2,2));
                
                funPha = @(x,k) funp(xx(x,:),kk(pos(k),:))-Us(x)*Vs(k)';
                ix = (1:N^2)'; ip = (1:numel(pos))';
                [Uph,Sph,Vph] = BF_lowrank(ix,ip,funPha,tol,5*mR,mR);
                U = Uph*Sph; V = Vph;
                U = [Us,U]; V = [Vs,V];
            end
            
            % low rank approximation of the smooth part
            funKsm = @(x,k) fun(xx(x,:),kk(pos(k),:))./(exp(1i*(U(x,1:rk)*V(k,1:rk)')));
            [Us,Ss,Vs] = BF_lowrank(ix,ip,funKsm,tol,5*mR,mR);
            Vs = Ss*Vs';
            
            % use NUFFT to evaluate the matvec for exp(1i*Us*Vs)
            % greengard's NUFFT
            if rk == 1 % use 1D NUFFT
                nufftfun = @(cj) nufft1d3(wlen,V(:,1),cj,1,1e-12,N^2,U(:,1));
                oscFun{cntw} = @(f) nufftfun(f);
            else if rk == 2 % use 2D NUFFT
                    nufftfun = @(cj) nufft2d3(wlen,V(:,1),V(:,2),cj,1,1e-12,N^2,U(:,1),U(:,2));
                    oscFun{cntw} = @(f) nufftfun(f);
                else if rk == 3 % use 3D NUFFT
                        nufftfun = @(cj) nufft3d3(wlen,V(:,1),V(:,2),V(:,3),cj,1,1e-12,N^2,U(:,1),U(:,2),U(:,3));
                        oscFun{cntw} = @(f) nufftfun(f);
                    end
                end
            end
            
            % apply the nonosc part
            rr = size(Us,2);
            maxRank3 = max(maxRank3,rr);
            appFun{cntw} = @(f) sum(Us.*oscFun{cntw}( (Vs.').*repmat(f(pos),[1,rr]) ),2);
        end
        FactorT3 = toc;
        
        tic;
        for cnt = 1:numTest
            yy = zeros(N^2,1);
            for cntw = 1:numW
                yy = yy + appFun{cntw}(f);
            end
        end
        ApplyT3 = toc/numTest;
        
        relerr3 = norm(yy(xs)-yext)/norm(yext);
        
        disp(['------------------------------------------']);
        disp(['N                 : ' num2str(N)]);
        disp(['maxRank           : ' num2str(maxRank3)]);
        disp(['Tolerance         : ' num2str(tol)]);
        disp(['Relative Error_2  : ' num2str(relerr3)]);
        disp(['Direct Time       : ' num2str(Td) ' s']);
        disp(['Factorization Time: ' num2str(FactorT3) ' s']);
        disp(['Applying Time     : ' num2str(ApplyT3) ' s']);
        disp(['------------------------------------------']);
    end
    
    %% old NUFFT approach
    numW = round(sqrt(N));
    oscFuno = cell(1,numW);
    appFuno = cell(1,numW);
    fundp1 = @(x,k)fun_phase_der1_2D(x,k);
    fundp2 = @(x,k)fun_phase_der2_2D(x,k);
    tic;
    Pos_angle(end/2+1,end/2+1) = 1e10;
    fc = fun(xx,[0,0]);
    pdc = sub2ind([N,N],N/2+1,N/2+1);
    maxRanko = 0;
    for cntw = 1:numW
        dw = 2*pi/numW; st = (cntw-1)*dw; ed = st + dw;
        idxMat = zeros(N,N);
        pos = find(Pos_angle<ed & Pos_angle>=st);
        idxMat(pos) = 1;
        pos = find(idxMat(:)>0);
        wlen = numel(pos);
        
        % low rank approximation of the phase
        agl = (st+ed)/2;
        kkc = [cos(agl),sin(agl)];
        der1 = fundp1(xx,kkc);
        der2 = fundp2(xx,kkc);
        U = [der1,der2];
        V = kk(pos,:);
        
        % low rank approximation of the smooth part
        ix = (1:N^2)'; ip = (1:numel(pos))';
        funKsm = @(x,k) fun(xx(x,:),kk(pos(k),:))./(exp(1i*(U(x,:)*V(k,:)')));
        [Us,Ss,Vs] = BF_lowrank(ix,ip,funKsm,tol,5*mR,mR);
        Vs = Ss*Vs';
        
        % use NUFFT to evaluate the matvec for exp(1i*Us*Vs)
        % greengard's NUFFT
        nufftfun = @(cj) nufft2d3(wlen,V(:,1),V(:,2),cj,1,1e-12,N^2,U(:,1),U(:,2));
        oscFuno{cntw} = @(f) nufftfun(f);
        
        % apply the nonosc part
        rr = size(Us,2);
        maxRanko = max(maxRanko,rr);
        appFuno{cntw} = @(f) sum(Us.*oscFuno{cntw}( (Vs.').*repmat(f(pos),[1,rr]) ),2);
    end
    FactorTo = toc;
    
    tic;
    for cnt = 1:numTest
        yy = zeros(N^2,1);
        for cntw = 1:numW
            yy = yy + appFuno{cntw}(f);
        end
        yy = yy + f(pdc)*fc;
    end
    ApplyTo = toc/numTest;
    
    relerro = norm(yy(xs)-yext)/norm(yext);
    
    disp(['------------------------------------------']);
    disp(['N                 : ' num2str(N)]);
    disp(['maxRank           : ' num2str(maxRanko)]);
    disp(['Tolerance         : ' num2str(tol)]);
    disp(['Relative Error_2  : ' num2str(relerro)]);
    disp(['Direct Time       : ' num2str(Td) ' s']);
    disp(['Factorization Time: ' num2str(FactorTo) ' s']);
    disp(['Applying Time     : ' num2str(ApplyTo) ' s']);
    disp(['------------------------------------------']);
    
    
    fprintf(fid,'%7d & %d & %.3e & %.2e & %.2e &  %d &  %.2e& %.2e& %.2e \\\\\n',...
        N,maxRanko,FactorTo,ApplyTo,relerro,maxRank,FactorT,ApplyT,relerr);
    fprintf(fid2,'%7d & %d & %.3e & %.2e & %.2e &  %d &  %.2e& %.2e& %.2e \\\\\n',...
        N,maxRanko,FactorTo,ApplyTo,relerro,maxRank3,FactorT3,ApplyT3,relerr3);
    
    
    
end
fprintf(fid,'\\ toprule \n');
fclose(fid);
fprintf(fid2,'\\ toprule \n');
fclose(fid2);





