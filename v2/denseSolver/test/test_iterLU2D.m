function test_iterLU2D(num)

% This is the test code for approximate LU factorization as a
% preconditioner in 2D for EFIE. We try to use the lower and upper
% triangular matrix of Z to approximate the LU factors of Z.
%
% Reference: A Hierarchical Butterfly LU Preconditioner for Two-Dimensional
% Electromagnetic Scattering Problems Involving Open Surfaces, preprint,
% 2019
% Copyright reserved.
% By Haizhao Yang, 2018, National University of Singapore

%
close all
vd = 10:15;
iterPre = zeros(size(vd));
iterNonPre = zeros(size(vd));
timeHSSBFapply = zeros(size(vd));
timeLUBF = zeros(size(vd));
timeLUBFapply = zeros(size(vd));
errIter = zeros(size(vd));
errDir = zeros(size(vd));
errRes = zeros(size(vd));
errResNon = zeros(size(vd));
errIterNon = zeros(size(vd));
method = 3; mr = 40;rk =8; tol = 1e-6;
restart = 50; tolsol = 1e-3; maxit = 100;
isShow = 0;
status = 2;
fname = ['./results/HBFLU/','mlu_',num2str(method),'_isRk_',num2str(rk),'_mr_',num2str(mr),'_isTol_',num2str(log10(1/tol)),'_',num2str(num)];
if isShow
    load([fname,'.mat']);
else
    for cnt = 1:numel(vd)
        N = 2^vd(cnt)+1;
        [N,num,rk]
        generate_Z;
        Afun = @(i,j) Z0fun(i,j)/nZ;
        
        ss = (1:N)'; tt = (1:N)';
        
        nps = max(3,ceil(log2(rk)));
        lsz = (2^nps)^2;
        
        % create the target linear system in a format of a function
        % handle for a fast matvec
        tic;
        [Afac,ZL,ZU] = HSSBF_RS_fwd(Afun,ss,tt,rk,tol,lsz,method);
        timeLUBF(cnt) = toc
        
        xt = f;
        Zfun = @(f) HSSBF_apply(Afac,f);
        tic;
        b = Zfun(xt);
        timeHSSBFapply(cnt) = toc
        
        Mfun = [];
        % solve Z*x=b without a preconditioner
        [x,flag,relres,iter,resvec] = gmres(Zfun,b,restart,tolsol,maxit,Mfun);
        iterNonPre(cnt) = (iter(1)-1)*restart+iter(2)
        errResNon(cnt) = relres
        errIterNon(cnt) = norm(x-xt)/norm(xt)
        
        % construct a preconditioner via triangular solvers
        Mfun = @(f )LUBF_sol(ZU,LUBF_sol(ZL, f,'L'),'U');
        tic;
        xdir = Mfun(b);
        timeLUBFapply(cnt) = toc
        errDir(cnt) = norm(xdir-xt)/norm(xt)
        % solve Z*x=b with the preconditioner Mfun
        [x4,flag4,relres4,iter4,resvec4] = gmres(Zfun,b,restart,tolsol,maxit,Mfun);
        iterPre(cnt) = (iter4(1)-1)*restart+iter4(2)
        
        errIter(cnt) = norm(x4-xt)/norm(xt)
        errRes(cnt) = relres4
       % save([fname,'.mat'],'vd','errIterNon','errResNon','errRes','iterPre','iterNonPre','timeHSSBFapply','timeLUBF','timeLUBFapply','errIter','errDir');
        
    end
end

% if 1
%     vs = 1:length(vd);
%     close all;
%     
%     pic = figure;
%     hold on;
%     h(1) = plot(vd(vs),log2(timeLUBF(vs)),'-^r','LineWidth',2);
%     h(2) = plot(vd(vs),log2(timeLUBFapply(vs)),'-*b','LineWidth',2);
%     vec = log2(2.^vd(vs).*vd(vs).*vd(vs));
%     h(3) = plot(vd(vs),log2(timeLUBF(vs(1)))+vec-vec(1),'--m','LineWidth',2);
%     vec = log2((2.^vd(vs)).^1.5);
%     h(4) = plot(vd(vs),log2(timeLUBF(vs(1)))+vec-vec(1),'--k','LineWidth',2);
%     vec = log2(2.^vd(vs).*vd(vs).*vd(vs));
%     h(5) = plot(vd(vs),log2(timeLUBFapply(vs(1)))+vec-vec(1),'--g','LineWidth',2);
%     vec = log2((2.^vd(vs)).^1.5);
%     h(6) = plot(vd(vs),log2(timeLUBFapply(vs(1)))+vec-vec(1),'--y','LineWidth',2);
%     legend(h,['HBFLU fac'],['HBFLU app'],'N log^2(N)','N^{1.5}','N log^2(N)','N^{1.5}','Location','EastOutside');
%     axis square;
%     xlabel('log_2(N)'); ylabel('log_2(time)');
%     set(gca, 'FontSize', 16);
%     b=get(gca);
%     set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
%     saveas(pic,[fname,'_tHBFLU.eps'],'epsc');
%     
%     pic = figure;
%     hold on;
%     h(1) = plot(vd(vs),iterNonPre(vs),'-^r','LineWidth',2);
%     h(2) = plot(vd(vs),iterPre(vs),'-*b','LineWidth',2);
%     legend(h,['no preconditioner'],['preconditioned'],'Location','NorthWest');
%     axis square;
%     xlabel('log_2(N)'); ylabel('iterations');
%     set(gca, 'FontSize', 16);
%     b=get(gca);
%     set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
%     saveas(pic,[fname,'_iterZ.eps'],'epsc');
% end
