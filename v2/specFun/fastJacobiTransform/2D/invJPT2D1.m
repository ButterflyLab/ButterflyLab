function [fun,rank1,rank2] = invJPT2D1(nts,da,db,tR,mR,tol,opt,R_or_N)
%% INPUTS:
%
%  nts    - the number of uniform smaples ts(bulit-in)
%  da,db  - polynomial parameters
%  tR     - p*mR, where p>5 sould be a oversampling parameter
%  mR     - maximum rank 
%  tol    - acquird accuracy
%  opt    - an option decides which lowrank approximation would be used (opt = 1 by default)
%           if opt >= 1, use Randomized sampling SVD [2]
%           if 0 <= opt < 1, use lowrank approximation based on applying chebyshev grids interpolation to samples [1], in this case, 
%              R_or_N should be smaller than 0.                   
%           if opt < 0, use lowrank approximation based on applying chebyshev
%           grids interpolation to samples and degrees [1], in this case, 
%              R_or_N should be smaller than 0.
%  R_or_N - a flag decides which method would be used (R_or_N = 1 by default)
%           if R_or_N > 0, the modified method [2] would be used
%           if R_or_N < 0, the original phase function method [1] would be used
%
%% OUTPUTS:
%
%  fun    -  a function handle computing the 2D inverse uniform Jacobi polynomial transform such that
%                  kron(J,J)*fun(c) = c(:), 
%            where  J(j,k) = M^(da,db)_(ts(j),k-1)*cos(psi^(da,db)_(ts(j),k-1)), 1=<j,k<=nts;
%  rank1  -  the rank of the lowrank approximation 
%  rank2  -  the rank of the lowrank approximation 
%
%%  For more details, please refer to 
%      
%      [1]. James Bremer and Haizhao Yang. Fast algorithms for Jacobi expansions via nonoscillatory
%      phase functions. arXiv:1803.03889 [math.NA], 2018.
%
%      [2]. James Bremer, Qiyuan Pang, Haizhao Yang. Fast Algorithms for the
%      Multi-dimensional Jacobi Polynomial Transform. arXiv:1901.07275 [math.NA], 2019.
%
%%  Copyright reserved by Qiyuan Pang, 25/1/2019

if nargin < 6
    disp('ERROR: ARGUMENTS DEFICIENT!')
    return
end

if  nargin == 6
    opt = 1;
    R_or_N = 1;
elseif nargin == 7
    R_or_N = 1;
end

if nargin > 8
    disp('ERROR: TOO MANY ARGUMENTS!')
    return
end

if R_or_N == 0
    disp('ERROR: ARGUMENT R_or_N SHOULD BE BIGGER THAN 0 OR SMALLER THAN 0!')
    return
end

if (opt<1)&&(R_or_N>0)
    disp('ERROR: WHEN ARGUMENT opt < 1, ARGUMENT R_or_N SHOULD BE SMALLER THAN 0!')
    return
end


   if nts < 2^12
       it = 10;
    else
       it = 28;
    end
nt = zeros(nts,1);
[ts,wghts] = getts(nt,da,db);
nu = [it:nts-1]';
xs = mod(floor(ts*nts/2/pi),nts)+1;

ts1 = ts;
ts2 = ts;
Y = [1:nts]';
X = zeros(nts,1);
S = ones(nts,1);
for i = 1:nts
    X(i) = xs(i);   
end
P1 = sparse(X,Y,S,nts,nts,nts);
P2 = P1;
% PP = kron(P2,P1);

if opt >= 1
    JTM = @(ts,nu)interpjac1(nt,ts,nu,da,db,R_or_N);
    [U,V] = lowrank(nts,JTM,ts,nu,tol,tR,mR);
   
    %U = repmat(sqrt(wghts),1,size(U,2)).*U;
    V = conj(V);
elseif 0 <= opt && opt<1
    [U,V] = ID_Cheby(nts,ts,nu,da,db,tol,1,R_or_N,tR,mR);
    U = repmat(sqrt(wghts),1,size(U,2)).*U;
    %U = diag(sqrt(wghts))*U;
elseif opt < 0
    [U,V] = ID_Cheby(nts,ts,nu,da,db,tol,-1,R_or_N,tR,mR);
    U = repmat(sqrt(wghts),1,size(U,2)).*U;
    %U = diag(sqrt(wghts))*U;
end
U1 = U;
U2 = U;
V1 = V;
V2 = V;
rank1 = size(U1,2);
rank2 = size(U2,2);
V1 = [zeros(it,rank1);V1];
V2 = [zeros(it,rank2);V2];
% UU = kron(U1,U2);
% VV = kron(V1,V2);
rank = rank1*rank2;

%vals1 = diag(sqrt(wghts))*jacrecur(nts,ts,it-1,da,db);
vals1 = jacrecur(nts,ts,wghts,it-1,da,db);
vals1 = vals1.';
%vals2 = diag(sqrt(wghts))*jacrecur(nts,ts,it-1,da,db);
vals2 = jacrecur(nts,ts,wghts,it-1,da,db);
vals2 = vals2.';
% vals = kron(vals1,vals2);

UU2 = reshape(repmat(U2,nts,1),nts,nts*rank2);
VV2 = reshape(repmat(V2,nts,1),nts,nts*rank2);
UU1 = reshape(repmat(U1,nts,1),nts,nts*rank1);
VV1 = reshape(repmat(V1,nts,1),nts,nts*rank1);

if  R_or_N > 0
    fun = @(c)invJacPT2d1(c);
else
    ex = exp(1i*nts/2*ts);
    U = U.*repmat(ex,1,rank);
    fun = @(c)invJacPT2d2(c);
end

    function y = invJacPT2d1(c)

        c = repmat(wghts,1,nts).*c.*repmat(wghts.',nts,1);
        y = zeros(nts,nts);
        y(1:it,:) = vals2*c;

        
        cc = repmat(c,1,rank2);
        dd = UU2.*cc;
%         dd = [];
%         for i = 1:rank2
%             dd = [dd d((i-1)*nts+1:i*nts,:)];
%         end
        dd = P2*dd;
        fftc = conj(fft(conj(dd)));
        
%         vv2 = [];
%         for i = 1:rank2
%             vv2 = [vv2 VV2((i-1)*nts+1:i*nts,:)];
%         end
        y2 = VV2.*fftc;
        y2 = reshape(y2,nts^2,rank2);
        y2 = reshape(sum(y2,2),nts,nts);
        y = y + real(y2);
        
        c = y.';
        
        y = zeros(nts,nts);
        y(1:it,:) = vals1*c;

        
        cc = repmat(c,1,rank1);
        dd = UU1.*cc;
%         dd = [];
%         for i = 1:rank1
%             dd = [dd d((i-1)*nts+1:i*nts,:)];
%         end
        dd = P1*dd;
        fftc = conj(fft(conj(dd)));
%         vv1 = [];
%         for i = 1:rank1
%             vv1 = [vv1 VV1((i-1)*nts+1:i*nts,:)];
%         end
        y1 = VV1.*fftc;
        y1 = reshape(y1,nts^2,rank1);
        y1 = reshape(sum(y1,2),nts,nts);
        y = y + real(y1);
        y = y.';
        y = y(:);

    end

end
