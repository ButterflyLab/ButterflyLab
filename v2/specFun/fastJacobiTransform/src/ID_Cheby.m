function [U,V] = ID_Cheby(n,x,k,da,db,tol,opt,R_or_N,tR,mR)
% Compute decomposition A = U*V.' via ID approximation A(:,rd) ~ A(:,sk)*T. 
% A =fun(rs,cs)
% x -- sample
% k -- degree
%
% tol -- accuracy
% da,db -- parameters of jacobi polynomial
% The precision is specified by tol and the rank is given by 'mR'; 
% 
% opt - whether use adaptive cheb interpolation for both sample dimension and degree dimension
%       if opt = 0, just use cheb interpolation for samples
%       if opt = 1, cheb interpolation for both samples and degrees
%        
% R_or_N -- decides which form of the matrix to be approximated
%       if R_or_N > 0, use the form : M^(da,db)(v,t)exp(i(psi^(a,b)(v,t)-2pi/n*[t*n/2pi]v))
%       if R_or_N <= 0, use the form : M^(da,db)(v,t)exp(i(psi^(a,b)(v,t)-tv))
%
% mR -- maximum rank
% tR -- p*mR, where p should be a oversampling paramter
% Copyright 2018 Haizhao Yang, Qiyuan Pang


%dd      = n;
%dd      = min(0.10,1/dd);
%dd      = log2(dd);
%nints   = ceil(-dd)+1;
%nints   = 2*nints;
%5kktol = ceil((log2(sqrt(pi/2)/tol)-0)/(4-log2(pi*exp(1))))*ceil(log2(n));
%kk = max(20*(ceil(log2(n))+1),kktol);


   k1 = 16;
   chebdata1 = cos([k1-1:-1:0]*pi/(k1-1))';

   dd         = 1/n;
   dd         = min(0.01,dd);

   dd         = log(dd)/log(2);
   nints      = ceil(-dd)+1;
   nintsab    = 2*nints;

   nints0  = nintsab/2;
   dd      = 2.0;
   
   ab = zeros(2,nintsab);
   for int=1:nints0
      ab(1,int) = dd^(-nints0+int-1);
      ab(2,int) = dd^(-nints0+int);
   end 

   ab = pi/2*ab;

   for int=1:nints0
      ab(1,nintsab-int+1) = pi - ab(2,int);
      ab(2,nintsab-int+1) = pi - ab(1,int);
   end 
   idx = 0;
   ts = zeros(nintsab*k1,1);
   for intab=1:nintsab
       a = ab(1,intab);
       b = ab(2,intab);
       for i=1:k1
           t = chebdata1(i)*(b-a)/2 + (b+a)/2;
           idx = idx+1;
           ts(idx)    = t;
       end
   end
   
   
   if  opt > 0
       nu = k;
   else
       k2 = 24;
       chebdata2 = cos([k2-1:-1:0]*pi/(k2-1))';
       cd0 = zeros(2,1000);
       nintscd = 0;
       if n < 2^12
          nintscd        =  nintscd+1;
          cd0(1,nintscd) =  9.0;
          cd0(2,nintscd) =  27.0;
       else
          nintscd        =  nintscd+1;
          cd0(1,nintscd) =  27.0;
          cd0(2,nintscd) =  81.0;
       end


       while cd0(2,nintscd) < n
             nintscd        = nintscd+1;
             cd0(1,nintscd) =  cd0(2,nintscd-1);
             cd0(2,nintscd) =  cd0(2,nintscd-1)*3.0;
       end

       cd0(2,nintscd) = min(n,cd0(2,nintscd));

       cd = cd0(:,1:nintscd); 
       dnus = zeros(nintscd*k2,1);
       for intcd=1:nintscd
           c = cd(1,intcd);
           d = cd(2,intcd);

           for i=1:k2
               dnu                = chebdata2(i)*(d-c)/2 + (d+c)/2;
               idx                = i+(intcd-1)*k2;
               dnus(idx)  = dnu;
           end
       end
       nu = dnus;
   end

%%%%%%%%%% consturct lowrank approximation for the matrix obtained via chebyshev grids
if opt > 0
    nt = zeros(n,1);
    A = interpjac1(nt,ts,nu,da,db,R_or_N);
    sigma = pi*chebdata1+pi;
    A = [A;exp(1i*sigma*nu'/n)];
    [~,R,E] = qr(A',0);
    rr = find( abs(diag(R)/R(1)) > tol, 1, 'last');
    rr = min(rr,mR);
    sk = E(1:rr);
    rd = E(rr+1:end);
    T = R(1:rr,1:rr)\R(1:rr,rr+1:end);
    T = T';
    V1 = A(sk,:);
    V1 = V1.';
    U1 = zeros(size(ts,1)+k1,rr);
    for i = 1:size(ts,1)+k1
        flag = find(sk == i);
        if ~isempty(flag)
            U1(i,flag) = 1;
        else
            flag1 = find(rd == i);
            U1(i,:) = T(flag1,:);
        end
    end
   
%%%%%%%%% construct right factor U in fun(x,k) = U*V.'
    w = [1/2 (-1).^[1:k1-2] 1/2*(-1)^(k1-1)]';
    S = barcycheby(x,chebdata1,w,ab);
    U = S*U1(1:size(ts,1),:);
    rx = x - floor(x*n/2/pi)*2*pi/n;
    [rx1,I] = sort(rx);
    J = zeros(size(rx,1),1);
    for i = 1:size(rx,1)
        J(I(i)) = i;
    end

    SS = barcycheby(rx1,chebdata1,w,[0.0;2*pi/n]);
    SS = SS(J,:);
    U = (SS*U1(size(ts,1)+1:end,:)).*U;
%%%%%%%%%% construct left factor V in fun(x,k) = U*V.'
    V = V1;
else
    nt = zeros(n,1);
    A = interpjac1(nt,ts,nu,da,db,R_or_N);
    sigma = pi*chebdata1+pi;
    B = exp(1i*sigma*nu'/n);
    A = [A;B];
    [~,R,E] = qr(A,0);
    rr = find( abs(diag(R)/R(1)) > tol, 1, 'last');
    rr = min(mR,rr);
    sk = E(1:rr);
    rd = E(rr+1:end);
    T = R(1:rr,1:rr)\R(1:rr,rr+1:end);
    U1 = A(:,sk);
    V1 = zeros(rr,size(nu,1));
    for i = 1:size(nu,1)
        flag = find(sk == i);
        if ~isempty(flag)
            V1(flag,i) = 1;
        else
            flag1 = find(rd == i);
            V1(:,i) = T(:,flag1);
        end
    end
    V1 = V1.';
    %[U1,S1,V1] = svdtrunc([A;B],mR,tol);
    %[U1,S1,V1] = svd([A;B]);
    %U1 = U1(:,1:rr);
    %V1 = V1(:,1:rr);
    %S1 = S1(1:rr,1:rr);
    w = [1/2 (-1).^[1:k1-2] 1/2*(-1)^(k1-1)]';
    S = barcycheby(x,chebdata1,w,ab);
    %U = S*U1(1:size(ts,1),:);
    U = S*U1(1:size(ts,1),:);
    rx = x - floor(x*n/2/pi)*2*pi/n;
    [rx1,I] = sort(rx);
    J = zeros(size(rx,1),1);
    for i = 1:size(rx,1)
	J(I(i)) = i;
    end
    
    SS = barcycheby(rx1,chebdata1,w,[0.0;2*pi/n]);
    SS = SS(J,:);

    %[~,R,E] = qr([A;B],0);
    %rr = find( abs(diag(R)/R(1)) > tol/10, 1, 'last');
    %rr = min(mR,rr);
    %rr = size(S1,1);
    %sk = E(1:rr);
    %rd = E(rr+1:end);
    %T = R(1:rr,1:rr)\R(1:rr,rr+1:end);
    %U2 = B(:,sk);
    
    U = (SS*U1(size(ts,1)+1:end,:)).*U;
    %U = (SS*U2).*U;
    w = [1/2 (-1).^[1:k2-2] 1/2*(-1)^(k2-1)]';
    P = barcycheby(k,chebdata2,w,cd);
    %V = P*conj(V1)*S1;
    %TE = interpjac1(nt,x,k,da,db,R_or_N);
    %norm(TE-U*V.')/norm(TE)
    V = P*V1;
end

end
