% Reference: Multidimensional Phase Recovery and Interpolative 
% Decomposition Butterfly Factorization, preprint, 2019.
% Copyright 2018 by Ze Chen, Juan Zhang, Ken Ho, and Haizhao Yang

clear; clc;
close all hidden;
rng('shuffle');

% -----input-----
datacase = '2D'; %'2D' '3D' 'Hel'
tol = 1e-9;
datatry = 10;
nn = 8;
tau = 1/4;
% ---------------

tstart = tic;

switch datacase
    case '2D'
        Nsq = [16^2, 16.^(2:6)];
        mR = 20;
        rk = 30;
    case '3D'
        Nsq = [8^4, 8.^(4:7)];
        mR = 3;
        rk = 80;
    case 'Hel'
        refinement = [4, 4:1:8];
        Nsq = 10*4.^(refinement-1);
        mR = 50;
        rk = 50;
end

startT = char(datetime('now'));
fid = fopen(['./results/' datacase '_' startT '.log'],'a+');
fprintf(fid, ['\n\n Starting from ', startT, '\n']);
fprintf(fid, [' case: ', datacase, '\tepsilon: ', num2str(tol), ...
    '\tr_epsilon: ', num2str(mR), '\tk_epsilon: ', num2str(rk), '\n']);
fprintf(fid, ['N & BF Error & Ker Error & Dis Row & Dis Col & '...
    'Ord Time & Rec Time & Fac Time & App Time & Dir/App Time\n']);

cm = length(Nsq);
LR = zeros(cm, datatry); Dis_r = LR; Dis_c = LR; Err_K = LR; Err_bf = LR;
T_order = LR; T_rec = LR; T_fac = LR; T_app = LR; T_dir = LR;

proc = 0;
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
                h = sqrt(N)/10;
                fun = @(xx,yy) fun_kernal_hel(h,xx,yy);
        end

        % Kernal function
        fn = @(i) fun(Sx,Sy(i,:));
        fnt = @(i) fun(Sx(i,:),Sy).';
        
        proc = ((c-1)*datatry+i-1)/datatry/cm;
        fprintf(['%2.2f%% \t N = %i \t %2.0f/%2.0f \t ',...
            char(datetime('now')),'\n'], proc*100, Nsq(c), i, datatry);
        
        % Ordering
        tord = tic;
        P1 = BF_RecoveryPath(Sx); P2 = BF_RecoveryPath(Sy);
        T_order(c,i) = toc(tord);

        % Recovery
        trec = tic;
        [U,S,V,Dr,Dc] = MIDBF_Lowrank_Phase(N,N,fn,fnt,tol,mR,mR,tau*2*pi,P1,P2);
        T_rec(c,i) = toc(trec);
        Dis_r(c,i) = Dr; Dis_c(c,i) = Dc;
        LR(c,i) = sum(diag(S)>tol*S(1));
        
        % Error
        NC = min(256, N);
        pos1 = BF_RandSample(N, NC);
        pos2 = BF_RandSample(N, NC);
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
        
        save(['./results/' datacase '_' startT '.mat'],'Nsq','Err_bf','Err_K',...
        'Dis_r','Dis_c','T_order','T_rec','T_fac','T_app','T_dir','LR');

        MIDBF_scaling([datacase '_' startT], 2, 2, 0)

    end
    
    fprintf(fid, ['%i & %.2e & %.2e & %.1f & %.1f & %.2e & %.2e & ', ...
        '%.2e & %.2e & %.2e \\\\\n'], Nsq(c), BF_trimdata(Err_bf(c,:)), ...
        BF_trimdata(Err_K(c,:)), mean(Dis_r(c,:)), mean(Dis_c(c,:)), ...
        BF_trimdata(T_order(c,:)), BF_trimdata(T_rec(c,:)), ...
        BF_trimdata(T_fac(c,:)), BF_trimdata(T_app(c,:)), ...
        BF_trimdata(T_dir(c,:)./T_app(c,:)));
end

T = toc(tstart) / 60;
fprintf(fid, ' Running time = %.2f min\n', T);
fprintf('Finished. Running time = %.2f min\n', T);
