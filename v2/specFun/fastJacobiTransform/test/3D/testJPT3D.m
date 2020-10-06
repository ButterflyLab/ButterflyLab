%% test JPT3D.m
MultiDimJacobi_startup
format long
num=20;
da=0.40;
db=0.40;
tol=1e-8
str1='size';
str2='our_rank1';
str3='our_rank2';
str4='our_time';
str5='our_rank3';
%str6='ratio_our/nyu';
str7='error_our';
%str8='error_nyu';
str9='dir_time';
str10='fac_time';
%str11='error_cheb';
fprintf('\n');
fprintf('start 3D uniform Jacobi polynomial transform test:');
fprintf('\n');
fprintf('da = %1.2f,db = %1.2f\n',da,db);
%fprintf('%-6s%-11s%-11s%-11s%-15s%-15s%-15s%-15s%-15s%-14s%-10s\n',str1,str10,str2,str3,str4,str5,str6,str11,str7,str8,str9);
fprintf('%-6s%-11s%-11s%-11s%-15s%-15s%-15s%-15s\n',str1,str2,str3,str5,str7,str4,str9,str10);
%funnyu = @(rs,cs,n,da,db,ts,nu)funnyu1d(rs,cs,n,da,db,ts,nu);
%funour = @(rs,cs,n,da,db,ts,nu)funour1d(rs,cs,n,da,db,ts,nu);
vd = [4.4:0.2:10];
es = length(vd);
rank1 = zeros(es,1);
rank2 = zeros(es,1);
rank3 = zeros(es,1);
errorour = zeros(es,1);
timeour = zeros(es,1);
timefac = zeros(es,1);
timedir = zeros(es,1);
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

    [ts,wghts] = getts(nt,da,db);
    xs = mod(floor(ts*nts/2/pi),nts)+1;
    ts1 = ts;
    ts2 = ts;
    ts3 = ts;
    rand('state',sum(100*clock));
    %ts1 = unique(rand(nts,1)*(pi-2/nts)+1/nts);
    %ts2 = unique(rand(nts,1)*(pi-2/nts)+1/nts);
    %ts3 = unique(rand(nts,1)*(pi-2/nts)+1/nts);
    nu = [it:nts-1]';
    n1 = randsample(nts*nts*nts,round(m));
    %n1 = [1:nts*nts]';
    %d = c;
    %tic;
    
    dirfun = directjac3d(nts,ts1,ts2,ts3,wghts,wghts,wghts,n1,da,db,c,'FULL');
    tic;
    result3 = dirfun(c);
    %vals = jacrecur(nts,ts,it-1,da,db);
    %J = interpjac1(nt,ts,nu,da,db,1);
    %J = [zeros(nts,it) J];
    %F = exp(1i*2*pi/nts*(xs-1)*[0:nts-1]);
    %J = J.*F;
    %J(:,1:it) = vals;
    %result4 = real(kron(kron(J,J),J)*c);
    timedir(ii) = toc;
    %norm(result3-result4(n1))/norm(result3)



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
 


%    [U1,V1]=lowrank(nts,funnyu,da,db,tol,tR,mR,ts,nu);
%    [U2,V2]=lowrank(nts,funour,da,db,tol,tR,mR,ts,nu);
%    rank1=size(U1,2);
%    rank2=size(U2,2);
%    ncol = size(c,2);

    tic
    for i = 1:num
    [fun,rank1(ii),rank2(ii),rank3(ii)] = JPT3D1(nts,da,db,tR,mR,tol,1,1);
    end
    timefac(ii)=toc/num;
    %P = U*V.';
    %P = [vals P];
    %norm(J-P)/norm(J)
    tic;
    for j=1:num
        result2 = fun(c);
    end
    timeour(ii)=toc/num;
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
    errorour(ii)=norm(result2-result3)/norm(result3);
    fprintf('\n  %-5d %-9d  %-9d  %-9d  %-1.6E   %-1.6E   %-1.6E   %-1.6E\n',nts,rank1(ii),rank2(ii),rank3(ii),errorour(ii),timeour(ii),timedir(ii),timefac(ii));
  
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
    vd = log2(round(2.^vd));
    ag = (log2(timeour(1))+log2(timefac(1))+log2(timedir(1)))/3;
    h(1) = plot(vd,4*vd-4*vd(1)+ag,'-r','LineWidth',2);
    h(2) = plot(vd,3*vd+2*log2(vd)-3*vd(1)-2*log2(vd(1))+ag,'--b','LineWidth',2);
    h(3) = plot(vd,3*vd+log2(vd)-3*vd(1)-log2(vd(1))+ag,'--k','LineWidth',2);
    h(4) = plot(vd,vd+log2(vd)-vd(1)-log2(vd(1))+ag,'--c','LineWidth',2);
    h(5) = plot(vd,vd-vd(1)+ag,'--m','LineWidth',2);
    h(6) = plot(vd,log2(timedir),'-*k','LineWidth',2);
    h(7) = plot(vd,log2(timeour),'-^g','LineWidth',2);
    h(8) = plot(vd,log2(timefac),'-^r','LineWidth',2);
    legend('n^4','n^3 log^2 n','n^3 log n','n log n','n','time dir','time app','time fac','Location','NorthWest');
    %title('3D uniform JPT, time');
    axis tight;
    xlabel('log_2(n)'); ylabel('log_{2}(time)');
    set(gca, 'FontSize', 20);
    b=get(gca);
    set(b.XLabel, 'FontSize', 20);set(b.YLabel, 'FontSize', 20);set(b.ZLabel, 'FontSize', 20);set(b.Title, 'FontSize', 20);
    saveas(pic,['testJPT3D_time.eps'],'epsc');
    hold off;
    pic1 = figure;
    hold on;
    h(1) = plot(vd,log10(errorour),'--k','LineWidth',2);
    legend('relerr','Location','NorthWest');
    %title('3D uniform JPT, relerr');
    axis tight;
    xlabel('log_2(n)'); ylabel('log_{10}(relerr)');
    set(gca, 'FontSize', 20);
    b=get(gca);
    set(b.XLabel, 'FontSize', 20);set(b.YLabel, 'FontSize', 20);set(b.ZLabel, 'FontSize', 20);set(b.Title, 'FontSize', 20);
    saveas(pic1,['testJPT3D_err.eps'],'epsc');
    hold off;
    pic2 = figure;
    hold on;
    h(1) = plot(vd,rank1.*rank2.*rank3,'--b','LineWidth',2);
    legend('rank','Location','NorthWest');
    %title('3D uniform JPT, rank');
    axis tight;
    xlabel('log_2(n)'); ylabel('rank');
    set(gca, 'FontSize', 20);
    b=get(gca);
    set(b.XLabel, 'FontSize', 20);set(b.YLabel, 'FontSize', 20);set(b.ZLabel, 'FontSize', 20);set(b.Title, 'FontSize', 20);
    saveas(pic2,['testJPT3D_rank.eps'],'epsc');
    hold off;
