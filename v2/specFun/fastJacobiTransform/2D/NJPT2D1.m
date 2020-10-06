function [fun,rank1,rank2] = NJPT2D1(nts,ts1,ts2,da,db,tR,mR,tol,opt,R_or_N)
%% INPUTS:
%
%  nts    - the number of the specific smaples ts1(ts2)
%  ts1    - the first specific samples given by users
%  ts2    - the second specific samples given by users
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
%  fun    -  a function handle computing the 2D forward nonuniform Jacobi polynomial transform such that
%                 fun(c(:)) = kron(J1,J2)*c(:), 
%            where  J1(j,k) = M^(da,db)_(ts1(j),k-1)*cos(psi^(da,db)_(ts1(j),k-1)), 1=<j,k<=nts;
%                   J2(j,k) = M^(da,db)_(ts2(j),k-1)*cos(psi^(da,db)_(ts2(j),k-1)), 1=<j,k<=nts;
%  rank1  -  the rank of the lowrank approximation of the matrix
%            corresponding to J1 
%  rank2  -  the rank of the lowrank approximation of the matrix
%            corresponding to J2
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

if nargin < 8
    disp('ERROR: ARGUMENTS DEFICIENT!')
    return
end

if  nargin == 8
    opt = 1;
    R_or_N = 1;
elseif nargin == 9
    R_or_N = 1;
end

if nargin > 10
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
%[ts,wghts] = getts(nt,da,db);
nu = [it:nts-1]';
%xs = mod(floor(ts*nts/2/pi),nts)+1;


if opt >= 1
    JTM1 = @(ts,nu)interpjac1(nt,ts,nu,da,db,R_or_N);
    [U1,V1] = lowrank(nts,JTM1,ts1,nu,tol,tR,mR);
    %U = diag(sqrt(wghts))*U;
    V1 = conj(V1);
    JTM2 = @(ts,nu)interpjac1(nt,ts,nu,da,db,R_or_N);
    [U2,V2] = lowrank(nts,JTM2,ts2,nu,tol,tR,mR);
    %U = diag(sqrt(wghts))*U;
    V2 = conj(V2);
elseif 0 <= opt && opt<1
    [U1,V1] = ID_Cheby(nts,ts1,nu,da,db,tol,1,R_or_N,tR,mR);
    [U2,V2] = ID_Cheby(nts,ts2,nu,da,db,tol,1,R_or_N,tR,mR);
    %U = diag(sqrt(wghts))*U;
elseif opt < 0
    [U1,V1] = ID_Cheby(nts,ts1,nu,da,db,tol,-1,R_or_N,tR,mR);
    [U2,V2] = ID_Cheby(nts,ts2,nu,da,db,tol,-1,R_or_N,tR,mR);
    %U = diag(sqrt(wghts))*U;
end
rank1 = size(U1,2);
rank2 = size(U2,2);
V1 = [zeros(it,rank1);V1];
V2 = [zeros(it,rank2);V2];
%UU = kron(U1,U2);
%VV = kron(V1,V2);
rank = rank1*rank2;

wghts1 = ones(nts,1);
wghts2 = ones(nts,1);
vals1 = jacrecur(nts,ts1,wghts1,it-1,da,db);
vals2 = jacrecur(nts,ts2,wghts2,it-1,da,db);
%vals = kron(vals1,vals2);

xs1 = mod(floor(ts1*nts/2/pi),nts)+1;
xs2 = mod(floor(ts2*nts/2/pi),nts)+1;
% xsub = zeros(nts*nts,1);
% for i = 1:nts
%     for j =1:nts
%         xsub((i-1)*nts+j) = (xs1(i)-1)*nts + xs2(j);
%     end
% end

VV2 = reshape(repmat(V2,nts,1),nts,nts*rank2);
UU2 = reshape(repmat(U2,nts,1),nts,nts*rank2);
VV1 = reshape(repmat(V1,nts,1),nts,nts*rank1);
UU1 = reshape(repmat(U1,nts,1),nts,nts*rank1);

if  R_or_N > 0
    fun = @(c)NJacPT2d1(c);
else
    ex = exp(1i*nts/2*ts);
    U = U.*repmat(ex,1,rank);
    fun = @(c)NJacPT2d2(c);
end

    function y = NJacPT2d1(c)
        y = vals2*c(1:it,:);

        
        cc = repmat(c,1,rank2);
        dd = VV2.*cc;

        fftc = nts*ifft(dd);
        fftc = fftc(xs2,:);
        
        y2 = UU2.*fftc;
        y2 = reshape(y2,nts^2,rank2);
        y2 = reshape(sum(y2,2),nts,nts);

        y = y + real(y2);

        c = y.';
        
        y = vals1*c(1:it,:);
        
        cc = repmat(c,1,rank1);
        dd = VV1.*cc;
        fftc = nts*ifft(dd);
        fftc = fftc(xs1,:);
        y1 = UU1.*fftc;
        y1 = reshape(y1,nts^2,rank1);
        y1 = reshape(sum(y1,2),nts,nts);
        y = y + real(y1);
	    y = y.';
        y = y(:);
    end

    function y = NJacPT2d2(c)
        y = zeros(nts,1);
        for i=1:rank
            cj = nufft1dIInyumex(ts,1,tol,V(:,i).*c);
            y = y + U(:,i).*cj;
        end
	y = real(y)./sqrt(wghts);
    end



end
