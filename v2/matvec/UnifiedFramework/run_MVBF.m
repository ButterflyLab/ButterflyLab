function run_MVBF(N, fioComp, fioCompt, NG, mR, tol, fid,type)
% NG - the number of interpolation grid points
% mR - the maximum rank for the phase and amplitude reconstruction
%
% Reference: H. Yang, A Unified Framework for Oscillatory Integral 
% Transform: When to use NUFFT or Butterfly factorization? preprint, 2018.

dim = 1;
tic;
% low-rank approximation of the phase function
[U,S,V] = MVBF_Lowrank_Phase(N,N,fioComp,fioCompt,tol,mR,mR,dim,pi/2,1);

% low-rank approximation of the amplitude function
[Ua,Sa,Va] = MVBF_Lowrank_Amp(N,N,fioComp,fioCompt,tol,mR,mR,dim);
recT = toc;

f = randn(N,1) + 1i*randn(N,1);

switch type
    case 1
        % IBF_multiscale_uniform
        V = S*V'; Va = Sa*Va';
        funk = @(t,s) exp(1i*U(t,:)*V(:,s));
        tt = 1:N; tt = tt(:);
        ss = 1:N; ss = ss(:);
        
        tic;
        [Factor,~] = IBF_multiscale_uniform(funk,tt,ss,NG,tol);
        FactorT = toc;
        
        rk = size(Ua,2);
        tic;
        yy = BF_multiscale_apply(Factor,(Va.').*repmat(f,[1,rk]));
        yy = sum(Ua.*yy,2);
        ApplyT = toc;
        RunT = FactorT + ApplyT;
    case 2 % CURBF
        V = S*V'; Va = Sa*Va';
        funk = @(t,s) exp(1i*U(t,:)*V(:,s));
        tt = 1:N; tt = tt(:);
        ss = 1:N; ss = ss(:);
        
        tic;
        [Factor,~] = CURBF(funk,tt,ss,NG,tol);
        FactorT = toc;
        
        rk = size(Ua,2);
        tic;
        yy = BF_apply(Factor,(Va.').*repmat(f,[1,rk]));
        yy = sum(Ua.*yy,2);
        ApplyT = toc;
        RunT = FactorT + ApplyT;
    case 3 % IBF_uniform with a simple partition dealing with the discontinuity of FIOs, but the gap grows with N, otherwise error is large
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
        
end
% check accuracy
ext = fioComp(f);
relerr = norm(ext-yy)/norm(ext);

fprintf(fid,'------------------------------------------\n');
fprintf(fid,'N                 : %4d\n', N);
fprintf(fid,'NG                : %4d\n', NG);
fprintf(fid,'mR                : %4d\n', mR);
fprintf(fid,'Tolerance         : %.3e\n', tol);
fprintf(fid,'BF Error          : %.3e\n', relerr);
fprintf(fid,'Running Time      : %.3e mins\n', RunT/60);
fprintf(fid,'Factorization Time: %.3e mins\n', FactorT/60);
fprintf(fid,'Reconstruct Time  : %.3e mins\n', recT/60);
fprintf(fid,'Applying Time     : %.3e s\n', ApplyT);
fprintf(fid,'------------------------------------------\n\n');

fprintf('BF Error          : %.3e\n', relerr);
end
