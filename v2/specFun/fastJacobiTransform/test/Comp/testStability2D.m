MultiDimJacobi_startup
format long
num=20;
da=0.75;
db=0.75;
tol=1e-14
str1='size';
str2='our_rank1';
str3='our_rank2';
str4='our_time';
%str5='nyu_time';
%str6='ratio_our/nyu';
str7='error_our';
%str8='error_nyu';
str9='dir_time';
str10='fac_time';
%str11='error_cheb';
fprintf('\n');
fprintf('start 2D uniform Jacobi polynomial transform stability test:');
fprintf('\n');
fprintf('da = %1.2f,db = %1.2f\n',da,db);
%fprintf('%-6s%-11s%-11s%-11s%-15s%-15s%-15s%-15s%-15s%-14s%-10s\n',str1,str10,str2,str3,str4,str5,str6,str11,str7,str8,str9);
fprintf('%-6s%-11s%-11s%-15s%-15s%-15s\n',str1,str2,str3,str7,str4,str10);
%funnyu = @(rs,cs,n,da,db,ts,nu)funnyu1d(rs,cs,n,da,db,ts,nu);
%funour = @(rs,cs,n,da,db,ts,nu)funour1d(rs,cs,n,da,db,ts,nu);
vd = [5:13];
es = length(vd);
rank1 = zeros(es,1);
rank2 = zeros(es,1);
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
    c = randn(nts,nts);
    d = c(:);

    [ts,wghts] = getts(nt,da,db);
    %ts = unique(rand(nts,1)*(pi-2/nts)+1/nts);
    xs = mod(floor(ts*nts/2/pi),nts)+1;
    nu = [it:nts-1]';
    %n1 = randsample(nts*nts,m);
    n1 = [1:nts*nts]';
    %d = c;
    tic;
    
%     result3 = directjac2d(nts,ts,ts,n1,da,db,d);
    %result1 = directinvjac2d(nts,ts,ts,wghts,wghts,n1,da,db,result3);
    %vals = jacrecur(nts,ts,it-1,da,db);
    %J = interpjac1(nt,ts,nu,da,db,1);
    %J = [zeros(nts,it) J];
    %F = exp(1i*2*pi/nts*(xs-1)*[0:nts-1]);
    %J = J.*F;
    %J(:,1:it) = vals;
    %result3 = real(kron(J,J)*c);
    timedir = nts*nts/m*(toc);
    %norm(result1-d)/norm(d)
    %norm(result4-result3(1:10))/norm(result3(1:10))



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
    [fun,rank1(ii),rank2(ii)] = JPT2D1(nts,da,db,tR,mR,tol,1,1);
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
    
    tic
    for i = 1:num
    [fun,rank1(ii),rank2(ii)] = invJPT2D1(nts,da,db,tR,mR,tol,1,1);
    end
    timefac(ii)=toc/num;
    %P = U*V.';
    %P = [vals P];
    %norm(J-P)/norm(J)
    tic;
    for j=1:num
        result = fun(reshape(result2,nts,nts));
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
    errorour(ii)=norm(result-c(:))/norm(c(:));
    fprintf('\n  %-5d %-9d  %-9d  %-1.6E   %-1.6E   %-1.6E\n',m,rank1(ii),rank2(ii),errorour(ii),timeour(ii),timefac(ii));
  
%    fprintf('\n   %-5d %-9d  %-9d  %-9d  %-1.6E   %-1.6E   %-1.6E   %-1.6E   %-1.6E   %-1.6E  %-1.6E\n',m,rank3,rank2,rank1,timeour,timenyu,timeratio,errorcheb,errorour,errornyu,timedir);
%    gc=imagesc(real(jacobi1(:,it+1:end)));
%    saveas(gc,'image13.jpg');
%    gf=imagesc(real(jacobi1(:,it+1:end)*1i));
%    saveas(gf,'image13i.jpg');
%    bf=imagesc(abs(jacobi1(:,it+1:end)));
%    saveas(bf,'image13a.jpg');
end
