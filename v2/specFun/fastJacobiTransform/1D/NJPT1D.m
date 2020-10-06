function [fun,rank] = NJPT1D(nts,ts,da,db,tR,mR,tol,opt,R_or_N)
%% INPUTS:
%
%  nts    - the number of the specific smaples ts(given by users)
%  ts     - the specific smaples given by users
%  da,db  - polynomial parameters
%  tR     - p*mR, where p>5 sould be a oversampling parameter
%  mR     - maximum rank 
%  tol    - acquird accuracy
%  opt    - an option decides which lowrank approximation would be used (opt = 0 by default)
%           if opt >= 1, use Randomized sampling SVD [2]
%           if 0 <= opt < 1, use lowrank approximation based on applying chebyshev grids interpolation to samples [1], in this case, 
%              R_or_N should be smaller than 0.                   
%           if opt < 0, use lowrank approximation based on applying chebyshev
%           grids interpolation to samples and degrees [1], in this case, 
%              R_or_N should be smaller than 0.
%  R_or_N - a flag decides which method would be used (R_or_N = -1 by default)
%           if R_or_N > 0, the modified method would [2] be used
%           if R_or_N < 0, the original phase function method [1] would be used
%
%% OUTPUTS:
%
%  fun    -  a function handle computing the 1D forward nonuniform Jacobi polynomial transform such that
%                 fun(c) = J*c, 
%            where  J(j,k) = M^(da,db)_(ts(j),k-1)*cos(psi^(da,db)_(ts(j),k-1)), 1=<j,k<=nts;
%  rank   -  the rank of the lowrank approximation 
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

if nargin < 7
    disp('ERROR: ARGUMENTS DEFICIENT!')
    return
end

if  nargin == 7
    opt = 0;
    R_or_N = -1;
elseif nargin == 8
    R_or_N = -1;
end

if nargin > 9
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
nu = [it:nts-1]';
xs = mod(floor(ts*nts/2/pi),nts)+1;
wghts = ones(nts,1);
nt = zeros(nts,1);

if  opt > 0
    JTM = @(ts,nu)interpjac1(nt,ts,nu,da,db,R_or_N);
    [U,V] = lowrank(nts,JTM,ts,nu,tol,tR,mR);
     V = conj(V);
elseif 0 <= opt && opt<1
    [U,V] = ID_Cheby(nts,ts,nu,da,db,tol,1,R_or_N,tR,mR);
elseif opt < 0
    [U,V] = ID_Cheby(nts,ts,nu,da,db,tol,-1,R_or_N,tR,mR);
end
rank = size(U,2);
V = [zeros(it,rank);V];

vals = jacrecur(nts,ts,wghts,it-1,da,db);

if  R_or_N > 0
    fun = @(c)NJacPT1d1(c);
else
    ex = exp(1i*nts/2*ts);
    U = U.*repmat(ex,1,rank);
    fun = @(c)NJacPT1d2(c);
end

    function y = NJacPT1d1(c)
        ncol = size(c,2);
        y = vals*c(1:it,:);

        VV = repmat(reshape(V,rank*nts,1),1,ncol);
        cc = repmat(c,rank,1);
        d = VV.*cc;
        dd = [];
        for i = 1:rank
            dd = [dd d((i-1)*nts+1:i*nts,:)];
        end
        fftc = nts*ifft(dd);
        fftc = fftc(xs,:);
        UU = repmat(reshape(U,rank*nts,1),1,ncol);
        uu = [];
        for i = 1:rank
            uu = [uu UU((i-1)*nts+1:i*nts,:)];
        end
        y2 = uu.*fftc;
        y2 = reshape(y2,nts*ncol,rank);
        y2 = reshape(sum(y2,2),nts,ncol);
        %y = y + real(y2);
        y = y + real(y2);
        
        
%         ncol = size(c,2);
%         d = repmat(V,1,ncol).*reshape(repmat(c,rank,1),nts,rank*ncol);
%         fftc = ifft(d);
%         fftc = fftc(xs,:);
%         y = nts*squeeze(sum(reshape(repmat(U,1,ncol).*fftc,nts,rank,ncol),2));
%         y = real(y);
%         y = y + vals*c(1:it-1,:);
    end
    
    function y = NJacPT1d2(c)
        y = zeros(nts,1);
        for i=1:rank
            cj = nufft1dIInyumex(ts,1,tol,V(:,i).*c);
            y = y + U(:,i).*cj;
        end
	    y = real(y);
        y = y + vals*c(1:it-1,:);
    end
end
