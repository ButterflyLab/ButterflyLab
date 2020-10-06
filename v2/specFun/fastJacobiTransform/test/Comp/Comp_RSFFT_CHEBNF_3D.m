%% This code compare the performance of the methods in [1] and [2] on computing 3D forward uniform (flag = 1) and nonuniform (flag = -1) tranforms.
%%  For more details, please refer to 
%      
%      [1]. James Bremer and Haizhao Yang. Fast algorithms for Jacobi expansions via nonoscillatory
%      phase functions. arXiv:1803.03889 [math.NA], 2018.
%
%      [2]. James Bremer, Qiyuan Pang, Haizhao Yang. Fast Algorithms for the
%      Multi-dimensional Jacobi Polynomial Transform. arXiv:1901.07275 [math.NA], 2019.
%
%%  Copyright reserved by Qiyuan Pang, 25/1/2019
MultiDimJacobi_startup
format long
flag = 1%don't change this value, now just works for flag > 0
num=10;
da=-0.80;
db=-0.80;
tol=1e-8
str1='size';
str2='RS_rank';
%str3='sCH_rank';
str6='sdCH_rank';
str4='RSapp_time';
%str5='sChapp_time';
str12='sdChapp_time';
str7='error_RS';
%str8='error_sCh';
str13='error_sdCh';
str9='dir_time';
str10='RSfac_time';
%str11='sChfac_time';
str14='sdChfac_time';
fprintf('\n');
fprintf('start RS SVD vs CHEB ID comparison 3D:');
fprintf('\n');
fprintf('da = %1.2f,db = %1.2f\n',da,db);
%fprintf('%-6s%-11s%-11s%-11s%-15s%-15s%-15s%-15s%-15s%-14s%-10s\n',str1,str10,str2,str3,str4,str5,str6,str11,str7,str8,str9);
fprintf('%-6s%-11s%-11s%-15s%-15s%-15s%-15s%-15s%-15s%-15s\n',str1,str2,str6,str4,str12,str7,str13,str9,str10,str14);
%funnyu = @(rs,cs,n,da,db,ts,nu)funnyu1d(rs,cs,n,da,db,ts,nu);
%funour = @(rs,cs,n,da,db,ts,nu)funour1d(rs,cs,n,da,db,ts,nu);
vd = [4.4:0.2:8];
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
    nts=round(2^m);
    if nts < 2^12
       it = 10;
    else
       it = 28;
    end

    nt=zeros(nts,1);
    c = randn(nts,nts,nts);
    d = c(:);

    if flag > 0
       [ts,wghts] = getts(nt,da,db);
       ts1 = ts;
       ts2 = ts;
       ts3 = ts;
       wghts1 = wghts;
       wghts2 = wghts;
       wghts3 = wghts;
    else
        %ts1 = unique(rand(nts,1)*(pi-2/nts)+1/nts);
        ts1 = abs(randn(nts,1));
        max1 = max(ts1);
        ts1 = unique(ts1/max1)*(pi-2/nts)+1/nts;
	wghts1 = ones(nts,1);
        %ts2 = unique(rand(nts,1)*(pi-2/nts)+1/nts);
        ts2 = abs(randn(nts,1));
        max2 = max(ts2);
        ts2 = unique(ts2/max2)*(pi-2/nts)+1/nts;
	wghts2 = ones(nts,1);
	%ts3 = unique(rand(nts,1)*(pi-2/nts)+1/nts);
	ts3 = abs(randn(nts,1));
        max3 = max(ts3);
        ts3 = unique(ts3/max3)*(pi-2/nts)+1/nts;
        wghts3 = ones(nts,1);
    end
    %ts = unique(rand(nts,1)*(pi-2/nts)+1/nts);
    nu = [0:nts-1]';
    n1 = randsample(nts*nts*nts,round(m));


    tic;

    result3 = directjac3d(nts,ts1,ts2,ts3,wghts1,wghts2,wghts3,n1,da,db,c);
    %size(v)
%    norm(result3)
    timedir = nts*nts*nts/m*(toc);




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
            [fun,r1,r2,r3] = JPT3D1(nts,da,db,tR,mR,tol,1,1);
        end
        rank1(ii) = r1*r2*r3;
        timefac1(ii)=toc/num;

        tic;
        for j=1:num
            result2 = fun(c);
        end
        timeour1(ii)=toc/num;

        errorour1(ii)=norm(result2(n1)-result3)/norm(result3);

        tic
        for i = 1:num
            [fun,r1,r2,r3] = JPT3D1(nts,da,db,tR,mR,tol,-1,-1);
        end
        rank2(ii) = r1*r2*r3;
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
            [fun,r1,r2,r3] = NJPT3D1(nts,da,db,tR,mR,tol,1,1);
        end
        rank1(ii) = r1*r2*r3;
        timefac1(ii)=toc/num;

        tic;
        for j=1:num
            result2 = fun(c);
        end
        timeour1(ii)=toc/num;

        errorour1(ii)=norm(result2(n1)-result3)/norm(result3);

        tic
        for i = 1:num
            [fun,r1,r2,r3] = NJPT3D1(nts,da,db,tR,mR,tol,-1,-1);
        end
        rank2(ii) = r1*r2*r3;
        timefac2(ii)=toc/num;

        tic;
        for j=1:num
            result2 = fun(c);
        end
        timeour2(ii)=toc/num;

        errorour2(ii)=norm(result2(n1)-result3)/norm(result3);

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

    fprintf('\n  %-5d %-9d %-9d %-1.6E   %-1.6E   %-1.6E   %-1.6E   %-1.6E   %-1.6E   %-1.6E\n',nts,rank1(ii),rank2(ii),timeour1(ii),timeour2(ii),errorour1(ii),errorour2(ii),timedir,timefac1(ii),timefac2(ii));

%    fprintf('\n   %-5d %-9d  %-9d  %-9d  %-1.6E   %-1.6E   %-1.6E   %-1.6E   %-1.6E   %-1.6E  %-1.6E\n',m,rank3,rank2,rank1,timeour,timenyu,timeratio,errorcheb,errorour,errornyu,timedir);
%    gc=imagesc(real(jacobi1(:,it+1:end)));
%    saveas(gc,'image13.jpg');
%    gf=imagesc(real(jacobi1(:,it+1:end)*1i));
%    saveas(gf,'image13i.jpg');
%    bf=imagesc(abs(jacobi1(:,it+1:end)));
%    saveas(bf,'image13a.jpg');
end
    figure('visible','off');
    vd = log2(round(2.^vd));
    pic = figure;
    hold on;
    %ag = (log2(timeour1(1))+log2(timeour2(1))+log2(timeour3(1)))/3;
    ag = (log2(timeour1(1))+log2(timeour2(1)))/2;
    h(1) = plot(vd,3*vd+log2(vd)-3*vd(1)-log2(vd(1))+ag,'--c','LineWidth',4);
    h(2) = plot(vd,3*vd+2*log2(vd)-3*vd(1)-2*log2(vd(1))+ag,'--k','LineWidth',4);
    h(3) = plot(vd,log2(timeour1),'-^r','LineWidth',2);
    %h(4) = plot(vd,log2(timeour2),'-^b','LineWidth',2);
    h(4) = plot(vd,log2(timeour2),'-^g','LineWidth',2);
    legend('n^3 log n','n^3 log^2 n','RS app','CHEB app','Location','bestoutside');
    %if flag > 0
    %   title('RS FFT vs CHEB NUFFT time, uni JPT');
    %else
    %    title('RS FFT vs CHEB NUFFT time, non JPT');
    %end
    axis tight;
    xlabel('log_2(n)'); ylabel('log_{2}(time)');
    set(gca, 'FontSize', 20);
    b=get(gca);
    set(b.XLabel, 'FontSize', 20);set(b.YLabel, 'FontSize', 20);set(b.ZLabel, 'FontSize', 20);set(b.Title, 'FontSize', 20);
    if flag > 0
       saveas(pic,['Comp_RSFFT_CHEBNF_uni1_3D080.eps'],'epsc');
    else
       saveas(pic,['Comp_RSFFT_CHEBNF_non1_3Di-080.eps'],'epsc');
    end
    hold off;
    pic = figure;
    hold on;
    %ag = (log2(timefac1(1))+log2(timefac1(1))+log2(timefac3(1)))/3;
    ag = (log2(timefac1(1))+log2(timefac2(1)))/2;
    h(1) = plot(vd,vd+2*log2(vd)-vd(1)-2*log2(vd(1))+ag,'--k','LineWidth',4);
    h(2) = plot(vd,vd+log2(vd)-vd(1)-log2(vd(1))+ag,'--c','LineWidth',4);
    h(3) = plot(vd,vd-vd(1)+ag,'--b','LineWidth',4);
    h(4) = plot(vd,log2(timefac1),'-xr','LineWidth',2);
    %h(4) = plot(vd,log2(timefac2),'-xb','LineWidth',2);
    h(5) = plot(vd,log2(timefac2),'-xg','LineWidth',2);
    legend('n log^2 n','n log n','n','RS fac','CHEB fac','Location','bestoutside');
    %if flag > 0
    %   title('RS FFT vs CHEB NUFFT time, uni JPT');
    %else
    %    title('RS FFT vs CHEB NUFFT time, non JPT');
    %end
    axis tight;
    xlabel('log_2(n)'); ylabel('log_{2}(time)');
    set(gca, 'FontSize', 20);
    b=get(gca);
    set(b.XLabel, 'FontSize', 20);set(b.YLabel, 'FontSize', 20);set(b.ZLabel, 'FontSize', 20);set(b.Title, 'FontSize', 20);
    if flag > 0
       saveas(pic,['Comp_RSFFT_CHEBNF_uni2_3D080.eps'],'epsc');
    else
       saveas(pic,['Comp_RSFFT_CHEBNF_non2_3D-080.eps'],'epsc');
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
    xlabel('log_2(n)'); ylabel('log_{10}(relerr)');
    set(gca, 'FontSize', 20);
    b=get(gca);
    set(b.XLabel, 'FontSize', 20);set(b.YLabel, 'FontSize', 20);set(b.ZLabel, 'FontSize', 20);set(b.Title, 'FontSize', 20);
    if flag > 0
       saveas(pic1,['CompFN_relerr_uni_3D080.eps'],'epsc');
    else
       saveas(pic1,['CompFN_relerr_non_3D-080.eps'],'epsc');
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
    xlabel('log_2(n)'); ylabel('rank');
    set(gca, 'FontSize', 20);
    b=get(gca);
    set(b.XLabel, 'FontSize', 20);set(b.YLabel, 'FontSize', 20);set(b.ZLabel, 'FontSize', 20);set(b.Title, 'FontSize', 20);
    if flag > 0
       saveas(pic2,['CompFN_rank_uni_3D080.eps'],'epsc');
    else
       saveas(pic2,['CompFN_rank_non_3D-080.eps'],'epsc');
    end

