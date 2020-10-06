%% This code compare the performance of the methods in [1] and [2] on computing 1D forward uniform (flag = 1) and nonuniform (flag = -1) tranforms.
%%  For more details, please refer to 
%      
%      [1]. James Bremer and Haizhao Yang. Fast algorithms for Jacobi expansions via nonoscillatory
%      phase functions. arXiv:1803.03889 [math.NA], 2018.
%
%      [2]. James Bremer, Qiyuan Pang, Haizhao Yang. Fast Algorithms for the
%      Multi-dimensional Jacobi Polynomial Transform. arXiv:1901.07275 [math.NA], 2019.
%
%%  Copyright reserved by Qiyuan Pang, 25/1/2019

format long
MultiDimJacobi_startup
flag = 1%don't change this value, now just works for flag > 0
num=10;
da=1.00;
db=1.00;
tol=1e-8
str1='size';
str2='RS_rank';
str3='CHEB_rank';
%str6='sdCH_rank';
str4='RSapp_time';
str5='CHEBapp_time';
%str12='sdChapp_time';
str7='error_RS';
str8='error_CHEB';
%str13='error_sdCh';
str9='dir_time';
str10='RSfac_time';
str11='CHEBfac_time';
%str14='sdChfac_time';
fprintf('\n');
fprintf('start RS SVD vs CHEB ID comparison:');
fprintf('\n');
fprintf('da = %1.2f,db = %1.2f\n',da,db);
%fprintf('%-6s%-11s%-11s%-11s%-15s%-15s%-15s%-15s%-15s%-14s%-10s\n',str1,str10,str2,str3,str4,str5,str6,str11,str7,str8,str9);
fprintf('%-6s%-11s%-11s%-15s%-15s%-15s%-15s%-15s%-15s%-15s\n',str1,str2,str3,str4,str5,str7,str8,str9,str10,str11);
%funnyu = @(rs,cs,n,da,db,ts,nu)funnyu1d(rs,cs,n,da,db,ts,nu);
%funour = @(rs,cs,n,da,db,ts,nu)funour1d(rs,cs,n,da,db,ts,nu);
vd = [8:24];
es = length(vd);
rank1 = zeros(es,1);
errorour1 = zeros(es,1);
timeour1 = zeros(es,1);
timefac1 = zeros(es,1);
rank2 = zeros(es,1);
errorour2 = zeros(es,1);
timeour2 = zeros(es,1);
timefac2 = zeros(es,1);
rank3 = zeros(es,1);
errorour3 = zeros(es,1);
timeour3 = zeros(es,1);
timefac3 = zeros(es,1);
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
    
    if flag > 0
       [ts,wghts] = getts(nt,da,db);
    else
        ts = unique(rand(nts,1)*(pi-2/nts)+1/nts);
        wghts = ones(nts,1);
    end
    %ts = unique(rand(nts,1)*(pi-2/nts)+1/nts);
    nu = [0:nts-1]';
    n1 = randsample(nts,m);

    d = c;
    tic;
    
    %[result3,t]=directjac1(nt,d,da,db,n1,ts,nu,wghts);
    %result3 = result3./sqrt(wghts(n1));
    %vals = jacrecur(nts,ts,it-1,da,db);
    %result3 = result3 + vals(n1,:)*d(1:it,:);
    %size(v)
%    norm(result3)    
    result3 =directjac1d(nts,ts,wghts,n1,da,db,c);
    timedir = nts/m*(toc);




    %xs=mod(floor(ts*nts/2/pi),nts)+1;
    %s=round(nts*ts);
    %gamma=norm(nts*ts-s,inf);
    %xi=log(log(10/tol)/gamma/7);
    %lw=xi-log(xi)+log(xi)/xi+0.5*log(xi)^2/xi^2-log(xi)/xi^2;
    %if m<10
    %   K=ceil(10*gamma*exp(lw));
    %elseif m<14
    %   K=ceil(12*gamma*exp(lw));
    %elseif m<18
    %   K=ceil(14*gamma*exp(lw));
    %elseif m<21
    %   K=ceil(16*gamma*exp(lw));
    %elseif m<24
    %   K=ceil(14*gamma*exp(lw));
    %elseif m<27
    %   K=ceil(15*gamma*exp(lw));
    %else
    %   K=ceil(17*gamma*exp(lw));
    %end
    %tR=K+2;
    %mR=K;
    p = 16;
    dd         = 1/nts;
    dd         = min(0.01,dd);

    dd         = log(dd)/log(2);
    nints      = ceil(-dd)+1;
    nints    = 2*nints;
    mR = ceil(2.0*log2(nts));
    tR = p*nints;


%    [U1,V1]=lowrank(nts,funnyu,da,db,tol,tR,mR,ts,nu);
%    [U2,V2]=lowrank(nts,funour,da,db,tol,tR,mR,ts,nu);
%    rank1=size(U1,2);
%    rank2=size(U2,2);
%    ncol = size(c,2);
    if  flag > 0
        tic
        for i = 1:num
            [fun,rank1(ii)] = JPT1D(nts,da,db,tR,mR,tol,1,1);
        end
        timefac1(ii)=toc/num;

        tic;
        for j=1:num
            result2 = fun(c);
        end
        timeour1(ii)=toc/num;
    
        errorour1(ii)=norm(result2(n1)-result3)/norm(result3);
    
        tic
        for i = 1:num
            [fun,rank2(ii)] = JPT1D(nts,da,db,tR,mR,tol,-1,-1);
        end
        timefac2(ii)=toc/num;

        tic;
        for j=1:num
            result2 = fun(c);
        end
        timeour2(ii)=toc/num;
    
        errorour2(ii)=norm(result2(n1)-result3)/norm(result3);
        
        %tic
        %for i = 1:num
        %    [fun,rank3(ii)] = JPT1D(nts,da,db,tR,mR,tol,-1,-1);
        %end
        %timefac3(ii)=toc/num;

        %tic;
        %for j=1:num
        %    result2 = fun(c);
        %end
        %timeour3(ii)=toc/num;
    
        %errorour3(ii)=norm(result2(n1)-result3)/norm(result3);
    else
        tic
        for i = 1:num
            [fun,rank1(ii)] = NJPT1D(nts,ts,da,db,tR,mR,tol,1,1);
        end
        timefac1(ii)=toc/num;

        tic;
        for j=1:num
            result2 = fun(c);
        end
        timeour1(ii)=toc/num;
    
        errorour1(ii)=norm(result2(n1)-result3)/norm(result3);
    
        tic
        for i = 1:num
            [fun,rank2(ii)] = NJPT1D(nts,ts,da,db,tR,mR,tol,0,-1);
        end
        timefac2(ii)=toc/num;

        tic;
        for j=1:num
            result2 = fun(c);
        end
        timeour2(ii)=toc/num;
    
        errorour2(ii)=norm(result2(n1)-result3)/norm(result3);
        
        tic
        for i = 1:num
            [fun,rank3(ii)] = NJPT1D(nts,ts,da,db,tR,mR,tol,-1,-1);
        end
        timefac3(ii)=toc/num;

        tic;
        for j=1:num
            result2 = fun(c);
        end
        timeour3(ii)=toc/num;
    
        errorour3(ii)=norm(result2(n1)-result3)/norm(result3);
    end
        
        
%    norm(result2)
%%%%%%%%%%%%%%%%%%%%%%%Greengard%%%%%%%%%%%%%%%%%
%    ex = exp(1i*nts/2*ts);
%    U1=U1.*repmat(ex,1,rank1);
%    tic;
%    for j=1:num
%        result1=zeros(nts,1);
%        for i=1:rank1
%            cj = nufft1dIInyumex(ts,1,tol,conj(V1(:,i)).*c);
%            result1 = result1 + U1(:,i).*cj;
%        end
%    end
%    timenyu=toc/num;
%    timeratio=timeour/timenyu;
%%    norm(result1)
%    [r,expvals,tss] = chebjacex(nt,da,db,tol);
%    rank3 = size(r,2);
%    xs=mod(floor(tss*nts/2/pi),nts)+1;
%    b = repmat(r,1,ncol).*reshape(repmat(c,rank3,1),nts,rank3*ncol);     
%    fftb = ifft(b);
%    fftb = fftb(xs,:);
%    result4 = nts*squeeze(sum(reshape(repmat(expvals,1,ncol).*fftb,nts,rank3,ncol),2));
%    errorcheb = norm(result4(n1)-result3)/norm(result3);
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         
    
    
%    error1=norm(result1-result2)/norm(result2)
%    errornyu=norm(result1(n1)-result3)/norm(result3);
    
    fprintf('\n  %-5d %-9d %-9d %-1.6E   %-1.6E   %-1.6E   %-1.6E   %-1.6E   %-1.6E   %-1.6E\n',m,rank1(ii),rank2(ii),timeour1(ii),timeour2(ii),errorour1(ii),errorour2(ii),timedir,timefac1(ii),timefac2(ii));
  
%    fprintf('\n   %-5d %-9d  %-9d  %-9d  %-1.6E   %-1.6E   %-1.6E   %-1.6E   %-1.6E   %-1.6E  %-1.6E\n',m,rank3,rank2,rank1,timeour,timenyu,timeratio,errorcheb,errorour,errornyu,timedir);
%    gc=imagesc(real(jacobi1(:,it+1:end)));
%    saveas(gc,'image13.jpg');
%    gf=imagesc(real(jacobi1(:,it+1:end)*1i));
%    saveas(gf,'image13i.jpg');
%    bf=imagesc(abs(jacobi1(:,it+1:end)));
%    saveas(bf,'image13a.jpg');
end
    figure('visible','off');
    pic = figure;
    hold on;
    %ag = (log2(timeour1(1))+log2(timeour2(1))+log2(timeour3(1)))/3;
    ag = (log2(timeour1(1))+log2(timeour3(1)))/2;
    h(1) = plot(vd,vd+log2(vd)-vd(1)-log2(vd(1))+ag,'--c','LineWidth',4);
    h(2) = plot(vd,vd+2*log2(vd)-vd(1)-2*log2(vd(1))+ag,'--k','LineWidth',4);
    h(3) = plot(vd,log2(timeour1),'-^r','LineWidth',2);
    %h(4) = plot(vd,log2(timeour2),'-^b','LineWidth',2);
    h(4) = plot(vd,log2(timeour2),'-^g','LineWidth',2);
    legend('N log N','N log^2 N','RS app','CHEB app','Location','bestoutside');
    %if flag > 0
    %   title('RS FFT vs CHEB NUFFT time, uni JPT');
    %else
    %    title('RS FFT vs CHEB NUFFT time, non JPT');
    %end
    axis tight;
    xlabel('log_2(N)'); ylabel('log_{2}(time)');
    set(gca, 'FontSize', 20);
    b=get(gca);
    set(b.XLabel, 'FontSize', 20);set(b.YLabel, 'FontSize', 20);set(b.ZLabel, 'FontSize', 20);set(b.Title, 'FontSize', 20);
    if flag > 0
       saveas(pic,['Comp_RSFFT_CHEBNF_uni1.eps'],'epsc');
    else 
       saveas(pic,['Comp_RSFFT_CHEBNF_non1.eps'],'epsc');
    end
    hold off;
    pic = figure;
    hold on;
    %ag = (log2(timefac1(1))+log2(timefac1(1))+log2(timefac3(1)))/3;
    ag = (log2(timefac1(1))+log2(timefac3(1)))/2;
    h(1) = plot(vd,vd+log2(vd)-vd(1)-log2(vd(1))+ag,'--c','LineWidth',4);
    h(2) = plot(vd,vd+2*log2(vd)-vd(1)-2*log2(vd(1))+ag,'--k','LineWidth',4);
    h(3) = plot(vd,log2(timefac1),'-xr','LineWidth',2);
    %h(4) = plot(vd,log2(timefac2),'-xb','LineWidth',2);
    h(4) = plot(vd,log2(timefac2),'-xg','LineWidth',2);
    legend('N log N','N log^2 N','RS fac','CHEB fac','Location','bestoutside');
    %if flag > 0
    %   title('RS FFT vs CHEB NUFFT time, uni JPT');
    %else
    %    title('RS FFT vs CHEB NUFFT time, non JPT');
    %end
    axis tight;
    xlabel('log_2(N)'); ylabel('log_{2}(time)');
    set(gca, 'FontSize', 20);
    b=get(gca);
    set(b.XLabel, 'FontSize', 20);set(b.YLabel, 'FontSize', 20);set(b.ZLabel, 'FontSize', 20);set(b.Title, 'FontSize', 20);
    if flag > 0
       saveas(pic,['Comp_RSFFT_CHEBNF_uni2.eps'],'epsc');
    else 
       saveas(pic,['Comp_RSFFT_CHEBNF_non2.eps'],'epsc');
    end
    hold off;
    pic1 = figure;
    hold on;
    h(1) = plot(vd,log10(errorour1),'-^r','LineWidth',2);
    %h(2) = plot(vd,log10(errorour2),'-^b','LineWidth',2);
    h(2) = plot(vd,log10(errorour2),'-^g','LineWidth',2);
    legend('RS relerr','CHEB relerr','Location','bestoutside');
    %if flag > 0
    %   title('relerr, uni JPT');
    %else
    %    title('relerr, non JPT');
    %end
    axis tight;
    xlabel('log_2(N)'); ylabel('log_{10}(relerr)');
    set(gca, 'FontSize', 20);
    b=get(gca);
    set(b.XLabel, 'FontSize', 20);set(b.YLabel, 'FontSize', 20);set(b.ZLabel, 'FontSize', 20);set(b.Title, 'FontSize', 20);
    if flag > 0
       saveas(pic1,['CompFN_relerr_uni.eps'],'epsc');
    else 
       saveas(pic1,['CompFN_relerr_non.eps'],'epsc');
    end
    hold off;
    pic2 = figure;
    hold on;
    h(1) = plot(vd,rank1,'-^r','LineWidth',2);
    %h(2) = plot(vd,rank2,'-^b','LineWidth',2);
    h(2) = plot(vd,rank2,'-^g','LineWidth',2);
    legend('RS rank','CHEB rank','Location','bestoutside');
    %if flag > 0
    %   title('rank, uni JPT');
    %else
    %    title('rank, non JPT');
    %end
    axis tight;
    xlabel('log_2(N)'); ylabel('rank');
    set(gca, 'FontSize', 20);
    b=get(gca);
    set(b.XLabel, 'FontSize', 20);set(b.YLabel, 'FontSize', 20);set(b.ZLabel, 'FontSize', 20);set(b.Title, 'FontSize', 20);
    if flag > 0
       saveas(pic2,['CompFN_rank_uni.eps'],'epsc');
    else 
       saveas(pic2,['CompFN_rank_non.eps'],'epsc');
    end
