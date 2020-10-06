%% test invJPT1D.m
MultiDimJacobi_startup
format long
num=20;
da=0.40;
db=0.40;
tol=1e-8
str1='size';
str2='our_rank';
str4='our_time';
str7='error_our';
str9='dir_time';
str10='time_fac';
fprintf('\n');
fprintf('start 1D inverse uniform Jacobi polynomial transform test:');
fprintf('\n');
fprintf('da = %1.2f,db = %1.2f\n',da,db);
fprintf('%-6s%-11s%-15s%-15s%-15s%-15s\n',str1,str2,str7,str4,str9,str10);
vd = [8:23];
es = length(vd);
rank = zeros(es,1);
errorour = zeros(es,1);
timeour = zeros(es,1);
timefac = zeros(es,1);
for ii=1:es
    m = vd(ii);
    nts=2^m;
    if nts < 2^12
       it = 10;
    else
       it = 28;
    end
    
    nt=zeros(nts,1);
    c = randn(nts,1);

    [ts,wghts] = getts(nt,da,db);
    %ts = unique(rand(nts,1)*(pi-2/nts)+1/nts);
    nu = [it:nts-1]';
    n1 = randsample(nts-it,m);
    d = c.*sqrt(wghts);
    
    tic;
    
    result3 =directinvjac1d(nts,ts,wghts,n1,da,db,c);
    %result3 = real(result3);
%    norm(result3)    
    timedir = nts/m*(toc);




    %xs=mod(floor(ts*nts/2/pi),nts)+1;
    %s=round(nts*ts);
    %gamma=norm(nts*ts-s,inf);
    %xi=log(log(10/tol)/gamma/7);
    %lw=xi-log(xi)+log(xi)/xi+0.5*log(xi)^2/xi^2-log(xi)/xi^2;
    %if m<10
    %   K=ceil(10*gamma*exp(lw));
    %elseif m<14
    %   K=ceil(11*gamma*exp(lw));
    %elseif m<18
    %   K=ceil(12*gamma*exp(lw));
    %elseif m<21
    %   K=ceil(13*gamma*exp(lw));
    %elseif m<24
    %   K=ceil(14*gamma*exp(lw));
    %elseif m<27
    %   K=ceil(15*gamma*exp(lw));
    %else
    %   K=ceil(17*gamma*exp(lw));
    %end
    p = 6;
    mR = ceil(1.5*log2(nts));
    tR = p*mR;

    tic
    for i = 1:num
    [fun,rank(ii)] = invJPT1D(nts,da,db,tR,mR,tol,1,1);
    end
    timefac(ii)=toc/num;

    tic;
    for j=1:num
        result2 = fun(c);
    end
    timeour(ii)=toc/num;

    errorour(ii)=norm(result2(n1)-result3)/norm(result3);
    fprintf('\n  %-5d %-9d  %-1.6E   %-1.6E   %-1.6E   %-1.6E\n',m,rank(ii),errorour(ii),timeour(ii),timedir,timefac(ii));
    
end
    figure('visible','off');
    pic = figure;
    hold on;
    ag = (log2(timeour(1))+log2(timefac(1)))/2;
    h(1) = plot(vd,vd+log2(vd)-vd(1)-log2(vd(1))+ag,'--k','LineWidth',2);
    h(2) = plot(vd,log2(timeour),'-^r','LineWidth',2);
    h(3) = plot(vd,log2(timefac),'-^g','LineWidth',2);
    legend('N log(N)','time app','time fac','Location','NorthWest');
    %title('1D inv uniform JPT, time');
    axis tight;
    xlabel('log_2(N)'); ylabel('log_{2}(time)');
    set(gca, 'FontSize', 20);
    b=get(gca);
    set(b.XLabel, 'FontSize', 20);set(b.YLabel, 'FontSize', 20);set(b.ZLabel, 'FontSize', 20);set(b.Title, 'FontSize', 20);
    saveas(pic,['testinvJPT1D_time.eps'],'epsc');
    hold off;
    pic1 = figure;
    hold on;
    h(1) = plot(vd,log10(errorour),'--k','LineWidth',2);
    legend('relerr','Location','NorthWest');
    %title('1D inv uniform JPT, relerr');
    axis tight;
    xlabel('log_2(N)'); ylabel('log_{10}(relerr)');
    set(gca, 'FontSize', 20);
    b=get(gca);
    set(b.XLabel, 'FontSize', 20);set(b.YLabel, 'FontSize', 20);set(b.ZLabel, 'FontSize', 20);set(b.Title, 'FontSize', 20);
    saveas(pic1,['testinvJPT1D_err.eps'],'epsc');
    hold off;
    pic2 = figure;
    hold on;
    h(1) = plot(vd,rank,'--b','LineWidth',2);
    legend('rank','Location','NorthWest');
    %title('1D inv uniform JPT, rank');
    axis tight;
    xlabel('log_2(N)'); ylabel('rank');
    set(gca, 'FontSize', 20);
    b=get(gca);
    set(b.XLabel, 'FontSize', 20);set(b.YLabel, 'FontSize', 20);set(b.ZLabel, 'FontSize', 20);set(b.Title, 'FontSize', 20);
    saveas(pic2,['testinvJPT1D_rank.eps'],'epsc');
