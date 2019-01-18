function scaling_IDBF(ex,method,rk,tol,vd,isShow,n0,rand_or_cheb,tt)
% ex - which example to test
% method - which N log N BF to use
% vd - set up problem sizes
% isShow - show results or compute results
% n0 - number of grid points per leaf
% tt - over sampling parameter in rand ID
% 
% Copyright 2018 by Haizhao Yang, National University of Singapore

% IDBF:
% rk = 30; tol = 1
% rk = 0; tol = 1e-6
% rk = 30; tol = 1e-6
% CURBF:
% rk = 6; tol = 1
% rk = 0; tol = 1e-6
% rk = 6; tol = 1e-6


w = warning ('on','all');
id = w.identifier;warning('off',id)
rmpath('folderthatisnotonpath')
%ex = 'B'; num = 1; vd = 8; isShow = 0;
if rk > 0
    isRk = 1; % use rk in the low-rank approximation
else
    isRk = 0;
end
if tol < 1
    isTol = 1; % use tol in the low-rank approximation
else
    isTol = 0;
end
fname = ['./results/IDBF/',ex,'_m_',num2str(method),'_fac_nn_',num2str(n0),'_isRk_',num2str(rk),'_isTol_',num2str(log10(1/tol)),'_',rand_or_cheb]; 
if nargin < 6, isShow = 0; end
if nargin < 5, vd = 8:17; end
if isShow
    load([fname,'.mat']);
    vs = 1:7;
else
    Thbfs = zeros(1,length(vd));
    Tbf = Thbfs;
    Tbfa = Thbfs;
    ebf = Thbfs;
    nzF = Thbfs;
    numTest = 5;
    for cnt = 1:length(vd)
        i = vd(cnt)
        
        switch ex
            case 'F'
                N = 2^i;
                ww = (-N/2:N/2-1)';
                rr = ((0:N-1)/N)';
                Z = @(x,k)(1+(sin(rr(x)/2)*cos(ww(k)/N/2)')/16).*fun_fio5_1D(rr(x),ww(k));
            case 'N'
                N = 2^i;
                ww = (-N/2:N/2-1)';
                ns = rand(N,1)/2;
                ww = ww + ns;
                rr = ((0:N-1)/N)';
                ns = rand(N,1)/N/2;
                rr = rr + ns;
                Z = @(x,k)funFT(rr(x),ww(k));
            case 'U'
                N = 2^i;
                ww = (-N/2:N/2-1)';
                rr = ((0:N-1)/N)';
                Z = @(x,k)funFT(rr(x),ww(k));
            case 'S'
                nv = 0; nu = 0;
                N = 2^i;
                w = ((1:N)+nv)*pi;
                ww = w(:);
                r = (0:N-1)/N;
                rr = r(:);
                Z = @(r,w) fun_Bessel(nu,rr(r),ww(w));
        end
        
        xx = (1:N)'; kk = (1:N)';
        NC = 256;
        xid = sort((1:NC)+round((N-NC)*rand(1,NC)));
        f = randn(N,1);
        yyt = Z(xx(xid),kk)*f;
        
        switch method
            case 1
                for ct = 1:numTest
                    tic;
                    Bfactor0 = CURBF2(Z,xx,kk,n0,rk,tol);
                    Tbf(cnt) = Tbf(cnt) + toc/numTest;
                    
                    tic;
                    yyl = BF_apply(Bfactor0,f);
                    Tbfa(cnt) = Tbfa(cnt) + toc/numTest;
                    relerr = norm(yyl(xid)-yyt)/norm(yyt);
                    ebf(cnt) = ebf(cnt) + relerr/numTest;
                    nzF(cnt) = nzF(cnt) + BF_nnz_show(Bfactor0)/numTest;
                end
            case 2
                for ct = 1:numTest
                    tic;
                    F = BF_IDBF(Z,xx,kk,n0,rk,tol,rand_or_cheb,tt);
                    Tbf(cnt) = Tbf(cnt) + toc/numTest;
                    tic;
                    yyl = BF_apply(F,f);
                    Tbfa(cnt) = Tbfa(cnt) + toc/numTest;
                    relerr = norm(yyl(xid)-yyt)/norm(yyt);
                    ebf(cnt) = ebf(cnt) + relerr/numTest;
                    nzF(cnt) = nzF(cnt) + BF_nnz_show(F)/numTest;
                end
        end
        ebf
        Tbf
        Tbfa
        nzF
        save([fname,'.mat'],'vd','Tbf','Tbfa','ebf','nzF');
    end
end

if 1
    vs = 2:length(vd);
    close all;
    
%     pic = figure;
%     hold on;
%     h(1) = plot(vd(vs),log2(Tbf(vs)),'-^r','LineWidth',2);
%     h(2) = plot(vd(vs),log2(Tbfa(vs)),'-*b','LineWidth',2);
%     vec = log2(2.^vd(vs).*vd(vs));
%     h(3) = plot(vd(vs),log2(100000*Tbf(vs(1)))+vec-vec(1),'--k','LineWidth',2);
%     legend(h,['BF fac'],['BF app'],'N log(N)','Location','NorthWest');
%     axis square;
%     xlabel('log_2(N)'); ylabel('log_2(time)');
%     set(gca, 'FontSize', 16);
%     b=get(gca);
%     set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
%     saveas(pic,[fname,'_rBF.eps'],'epsc');
    
    pic = figure;
    hold on;
    h(1) = plot(vd(vs),log2(Tbf(vs)),'-^r','LineWidth',2);
    h(2) = plot(vd(vs),log2(Tbfa(vs)),'-*b','LineWidth',2);
    vec = log2(2.^vd(vs).*vd(vs));
    h(3) = plot(vd(vs),log2(Tbf(vs(1)))+vec-vec(1),'--m','LineWidth',2);
    vec = log2(2.^vd(vs).*vd(vs).*vd(vs));
    h(4) = plot(vd(vs),log2(Tbf(vs(1)))+vec-vec(1),'--k','LineWidth',2);
    vec = log2(2.^vd(vs).*vd(vs));
    h(5) = plot(vd(vs),log2(Tbfa(vs(1)))+vec-vec(1),'--g','LineWidth',2);
    vec = log2(2.^vd(vs).*vd(vs).*vd(vs));
    h(6) = plot(vd(vs),log2(Tbfa(vs(1)))+vec-vec(1),'--y','LineWidth',2);
    legend(h,['BF fac'],['BF app'],'N log(N)','N log^2(N)','N log(N)','N log^2(N)','Location','EastOutside');
    axis square;
    xlabel('log_2(N)'); ylabel('log_2(time)');
    set(gca, 'FontSize', 16);
    b=get(gca);
    set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
    saveas(pic,[fname,'_tBF.eps'],'epsc');
    
    
    pic = figure;
    hold on;
    h(1) = plot(vd(vs),log10(ebf(vs)),'-^r','LineWidth',2);
    legend(h,['err BF'],'Location','NorthWest');
    axis square;
    xlabel('log_2(N)'); ylabel('log_{10}(error)');
    set(gca, 'FontSize', 16);
    b=get(gca);
    set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
    saveas(pic,[fname,'_err.eps'],'epsc');
    
    pic = figure;
    hold on;
    h(1) = plot(vd(vs),log2(nzF(vs)),'-^r','LineWidth',2);
    vec = log2(2.^vd(vs).*vd(vs));
    h(2) = plot(vd(vs),log2(nzF(vs(1)))+vec-vec(1),'--m','LineWidth',2);
    vec = log2(2.^vd(vs).*vd(vs).*vd(vs));
    h(3) = plot(vd(vs),log2(nzF(vs(1)))+vec-vec(1),'--k','LineWidth',2);
    legend(h,['nnz BF'],'N log(N)','N log^2(N)','Location','NorthWest');
    axis square;
    xlabel('log_2(N)'); ylabel('log_{2}(nnz)');
    set(gca, 'FontSize', 16);
    b=get(gca);
    set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
    saveas(pic,[fname,'_nnz.eps'],'epsc');
end



