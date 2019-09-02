% Reference: Multidimensional Phase Recovery and Interpolative 
% Decomposition Butterfly Factorization, preprint, 2019.
% Copyright 2018 by Ze Chen, Juan Zhang, Ken Ho, and Haizhao Yang

clear; clc;
close all hidden;

% -----input-----
datacase = 'Hel'; %'2D' '3D' 'Hel'
tol = 1e-9;
datatry = 10;
nn = 8;
% ---------------

tstart = tic;

switch datacase
    case '2D'
        Nsq = [64, 16.^(2:6)];
        mR = 10;
        rk = 30;
    case '3D'
        Nsq = [1024, 8.^(4:7)];
        mR = 5;
        rk = 80;
    case 'Hel'
        refinement = [3, 4:1:10];
        Nsq = 10*4.^(refinement-1);
        mR = 100;
        rk = 30;
end

fid = fopen(['./results/' datacase '.log'],'a+');
fprintf(fid,['\n\n Starting from ',char(datetime('now')),'\n']);
fprintf(fid,[' epsilon: ',num2str(tol),'\t', ...
    'r_epsilon: ',num2str(mR),'\t','k_epsilon: ',num2str(rk),'\n']);
fprintf(fid,['N & BF Error & Ker Error & ',...
    'Rec Time & Fac Time & App Time & Dir/App Time\n']);

cm = length(Nsq);
LR = zeros(cm,datatry); Err_K = LR; Err_bf = LR;
T_order = LR; T_rec = LR; T_fac = LR; T_app = LR; T_dir = LR;

proc = 0;
mp = 100/datatry+0.1;
fprintf('\n');

for c = 1 : cm
    
    for i = 1 : datatry
        
        switch datacase
            case '2D'
                N = Nsq(c);
                [X,Y] = meshgrid(unigrid(0,1/sqrt(N),1,'[)'));
                Sx = [X(:),Y(:)];
                Sy = Sx;
                fun = @(xx,yy) fun_fio_2D(xx,yy);
            case '3D'
                N = Nsq(c);
                Sx = rand(N,3);
                Sy = rand(N,3);
                fun = @(xx,yy) funFT(xx,yy);
            case 'Hel'
                part = 4; px = 4; py = 1;
                % Discretize sphere
                [~, sph, ~] = BF_Helmholtz(1, refinement(c));
                sph(:, 1) = sph(:, 1) / (2*pi) + 0.5;
                sph(:, 2) = sph(:, 2) / pi + 0.5;
                sph = sortrows(sph);
                Sx = sph; Sy = sph;
                N = length(sph);
                % Sample location
                x = round(N/part*px)-ceil(N/part-1)+1 : round(N/part*px);
                y = round(N/part*(py-1))+1 : ...
                    round(N/part*(py-1))+ceil(N/part-1);
                Sx = Sx(x, :); Sy = Sy(y, :);
                Nx = length(x); Np = length(y);
                N = min(Nx, Np);
                Nsq(c) = N;
                % Kernal function
                k = sqrt(N)/6;
                fun = @(xx,yy) fun_kernal_hel(k,xx,yy);
        end

        % Kernal function
        fn = @(i) fun(Sx,Sy(i,:));
        fnt = @(i) fun(Sx(i,:),Sy).';
        
        proc = ((c-1)*datatry+i-1)/datatry/cm;
        fprintf(['%2.2f%% \t N = %i \t %2.0f/%2.0f \t ',...
            char(datetime('now')),'\n'], ...
            proc*100, Nsq(c), i, datatry);
        
        % Ordering
        tord = tic;
        P1 = BF_RecoveryPath(Sx); P2 = BF_RecoveryPath(Sy);
        T_order(c,i) = toc(tord);

        % Recovery
        trec = tic;
        [U,S,V] = MIDBF_Lowrank_Phase(N,N,fn,fnt,tol,mR,mR,1/4*pi,P1,P2);
        T_rec(c,i) = toc(trec);
        LR(c,i) = sum(diag(S)>tol*S(1));
        
        % Error
        NC = min(256, N);
        pos1 = round(rand(1,NC)*(N-NC))+(1:NC);
        pos2 = round(rand(1,NC)*(N-NC))+(1:NC);
        K = fun(Sx(pos1,:), Sy(pos2,:));
        Err_K(c,i) = norm(exp(1i*U(pos1,:)*S*V(pos2,:)')-K)/norm(K);
        
        % optimal BF via random sampling and IBF_uniform
        f = randn(N,1) + 1i*randn(N,1);
        V = S*V';
        funk = @(x,y) exp(1i*U(x,:)*V(:,y));
        tic;
        Factor = BF_MIDBF_nonuniform(funk,Sx,Sy,nn,rk,tol,'cheb');
        T_fac(c,i) = toc;
        tic;
        yy = BF_apply(Factor, f);
        T_app(c,i) = toc;
        tic;
        Err_bf(c,i) = BF_check(N,fun,f,Sx,Sy,yy,NC);
        T_dir(c,i) = toc;
        T_dir(c,i) = T_dir(c,i) * N / NC;
        
    end
    
    fprintf(fid, '%i & %.2e & %.2e & %.2e & %.2e & %.2e & %.2e\n', ...
        Nsq(c),trimmean(Err_bf(c,:),mp),trimmean(Err_K(c,:),mp),...
        trimmean(T_order(c,:)+T_rec(c,:),mp),trimmean(T_fac(c,:),mp),...
        trimmean(T_app(c,:),mp),trimmean(T_dir(c,:)./T_app(c,:),mp));
end

save(['./results/' datacase '.mat'],'Nsq','Err_bf','Err_K','T_order',...
    'T_rec','T_fac','T_app','T_dir','LR','mp');

T = toc(tstart) / 60;
fprintf(fid, ' Running time = %.2f min\n', T);
fprintf('Finished. Running time = %.2f min\n', T);
