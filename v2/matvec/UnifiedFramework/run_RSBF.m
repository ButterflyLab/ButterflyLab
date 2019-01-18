function run_RSBF(N, func_name, NG, mR, tol, fid, method)
% NG - the number of interpolation grid points
% mR - the maximum rank for the phase and amplitude reconstruction
%
% Reference: H. Yang, A Unified Framework for Oscillatory Integral 
% Transform: When to use NUFFT or Butterfly factorization? preprint, 2018.

switch func_name
    case 'funF' % Fourier transform
        k = -N/2:(N/2-1);
        x = (0:N-1)/N;
        fun = @(x,k)funFT(x,k);
    case 'funIF' % Fourier transform
        k = -N/2:(N/2-1);
        x = (0:N-1)/N;
        fun = @(x,k)funIFT(x,k);
    case 'fun1' % Fourier integral operator
        k = -N/2:(N/2-1);
        x = (0:N-1)/N;
        fun = @(x,k)fun_fio5_1D(x,k);
    case 'fun2' % Fourier integral operator
        k = -N/2:(N/2-1);
        x = (0:N-1)/N;
        fun = @(x,k)fun_fio2_1D(x,k);
    case 'funH' % Hankel transformation
        k = 1:N;
        x = 1:N;
        fun = @(x,k)fun_Hankel(N,x,k);
end
kk = k(:);
xx = x(:);

switch method
    case 1
        % define sampling operators
        fn = @(y) fun(xx,kk(y));
        fnt = @(y) fun(xx(y),kk).';
        
        dim = 1;
        
        tic;
        % low-rank approximation of the phase function
        [U,S,V] = RSBF_Lowrank_Phase(N,N,fn,fnt,tol,mR,mR,dim,pi/2);
        
        % low-rank approximation of the amplitude function
        [Ua,Sa,Va] = RSBF_Lowrank_Amp(N,N,fn,fnt,tol,mR,mR,dim);
        recT = toc;
        
        %% optimal BF via random sampling and IBF_uniform
        f = randn(N,1) + 1i*randn(N,1);
        
        V = S*V'; Va = Sa*Va';
        funk = @(t,s) exp(1i*U(t,:)*V(:,s));
        tt = 1:N; tt = tt(:);
        ss = 1:N; ss = ss(:);
        
        tic;
        [Factor,~] = IBF_uniform(funk,tt,ss,NG,tol);
        FactorT = toc;
        
        rk = size(Ua,2);
        tic;
        yy = BF_apply(Factor,(Va.').*repmat(f,[1,rk]));
        yy = sum(Ua.*yy,2);
        ApplyT = toc;
        RunT = FactorT + ApplyT;
        
        % check accuracy
        NC = 256;
        pos1 = round(rand(1,NC)*(N-NC))+(1:NC);
        pos2 = round(rand(1,NC)*(N-NC))+(1:NC);
        K = fun(xx(pos1),kk(pos2));
        errAll = norm((Ua(pos1,:)*Va(:,pos2)).*exp(1i*U(pos1,:)*V(:,pos2))-K)/norm(K);
        errPh = norm(exp(1i*U(pos1,:)*V(:,pos2))-K./abs(K))/norm(K./abs(K));
        errAmp = norm((Ua(pos1,:)*Va(:,pos2))-abs(K))/norm(abs(K));
        
    case 2
        fn = @(x,y) fun(xx(x),kk(y));
        tt = 1:N; tt = tt(:);
        ss = 1:N; ss = ss(:);
        
        tic;
        [Factor,~] = IBF_uniform(fn,tt,ss,NG,tol);
        FactorT = toc;
        
        f = randn(N,1) + 1i*randn(N,1);
        tic;
        yy = BF_apply(Factor,f);
        ApplyT = toc;
        RunT = FactorT + ApplyT;
        NC = 256;
end

tic;
relerr = BF_check(N,fun,f,xx,kk,yy,NC);
Td = toc;
Td = Td*N/NC;

fprintf(fid,'------------------------------------------\n');
fprintf(fid,'N                 : %4d\n', N);
fprintf(fid,'NG                : %4d\n', NG);
fprintf(fid,'mR                : %4d\n', mR);
fprintf(fid,'Tolerance         : %.3e\n', tol);
fprintf(fid,'BF Error          : %.3e\n', relerr);
if method == 1
    fprintf(fid,'All Error         : %.3e\n', errAll);
    fprintf(fid,'Phase Error       : %.3e\n', errPh);
    fprintf(fid,'Amplitude Error   : %.3e\n', errAmp);
end
fprintf(fid,'Direct Time       : %.3e s\n', Td);
fprintf(fid,'Running Time      : %.3e mins\n', RunT/60);
fprintf(fid,'Factorization Time: %.3e mins\n', FactorT/60);
if method == 1
    fprintf(fid,'Reconstruct Time  : %.3e mins\n', recT/60);
end
fprintf(fid,'Applying Time     : %.3e s\n', ApplyT);
fprintf(fid,'------------------------------------------\n\n');

save(['./results/Factor_' func_name '_' num2str(N) '_' num2str(NG) '_' num2str(mR) '.mat'],'Factor','-v7.3');

end
