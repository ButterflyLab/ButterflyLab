% This code compare two N log N bottom-up butterfly factorization methods
% Run scaling_IDBF.m first to generate results. Then run this code to
% compare results.
% 
% Copyright 2018 by Haizhao Yang, National University of Singapore

w = warning ('on','all');
id = w.identifier;warning('off',id)
rmpath('folderthatisnotonpath')
exv = {'N','S','F'};
close all;
n0 = 8;
isRk = 0;
for cnt = 1:3
    clear Tbf Tbfa ebf;
    ex = exv{cnt};
    fname = ['./results/IDBF/',ex,'_comp_fac',num2str(n0)]; 
    method = 2;
    fname1 = ['./results/IDBF/',ex,'_m_',num2str(method),'_fac_nn_',num2str(n0),'_isRk_',num2str(isRk)]; 
    load([fname1,'.mat']);
    Tbf1 = Tbf; Tbfa1 = Tbfa; ebf1 = ebf;
    clear Tbf Tbfa ebf;
    method = 1;
    fname2 = ['./results/IDBF/',ex,'_m_',num2str(method),'_fac_nn_',num2str(n0),'_isRk_',num2str(isRk)]; 
    load([fname2,'.mat']);
    Tbf2 = Tbf; Tbfa2 = Tbfa; ebf2 = ebf;
    
    vd = 8:17;
    vs = 2:2:length(vd);
    
    pic = figure;
    hold on;
    h(1) = plot(vd(vs),log2(Tbf1(vs)),'^r','LineWidth',2);
    h(2) = plot(vd(vs),log2(Tbf2(vs)),'*b','LineWidth',2);
    vec = log2(2.^vd(vs).*vd(vs));
    h(3) = plot(vd(vs),log2(Tbf1(vs(1)))+vec-vec(1),'--k','LineWidth',2);
    h(4) = plot(vd(vs),log2(Tbfa1(vs)),'vm','LineWidth',2);
    h(5) = plot(vd(vs),log2(Tbfa2(vs)),'og','LineWidth',2);
    vec = log2(2.^vd(vs).*vd(vs));
    h(6) = plot(vd(vs),log2(Tbfa1(vs(1)))+vec-vec(1),'--k','LineWidth',2);
    legend(h,['IDBF fac'],['CURBF fac'],'N log(N)',['IDBF app'],['CURBF app'],'N log(N)','Location','NorthWest');
    axis square;
    xlabel('log_2(N)'); ylabel('log_2(time)');
    set(gca, 'FontSize', 16);
    b=get(gca);
    set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
    saveas(pic,[fname,'_time.eps'],'epsc');
    
    
    pic = figure;
    hold on;
    h(1) = plot(vd(vs),log10(ebf1(vs)),'-^r','LineWidth',2);
    h(2) = plot(vd(vs),log10(ebf2(vs)),'-*b','LineWidth',2);
    legend(h,['Err IDBF'],['Err CURBF'],'Location','NorthWest');
    axis square;
    xlabel('log_2(N)'); ylabel('log_{10}(error)');
    set(gca, 'FontSize', 16);
    b=get(gca);
    set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
    saveas(pic,[fname,'_err.eps'],'epsc');
end


