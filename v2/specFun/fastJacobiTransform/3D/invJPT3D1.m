function [fun,rank1,rank2,rank3] = invJPT3D1(nts,da,db,tR,mR,tol,opt,R_or_N)
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
%  fun    -  a function handle computing the 3D inverse uniform Jacobi polynomial transform such that
%                 kron(kron(J,J),J)*fun(c) = c(:), 
%            where  J(j,k) = M^(da,db)_(ts(j),k-1)*cos(psi^(da,db)_(ts(j),k-1)), 1=<j,k<=nts;
%  rank1  -  the rank of the lowrank approximation 
%  rank2  -  the rank of the lowrank approximation 
%  rank3  -  the rank of the lowrank approximation 
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
ts1 = ts;
ts2 = ts;
ts3 = ts;
%xs = mod(floor(ts*nts/2/pi),nts)+1;


if opt >= 1
    JTM1 = @(ts,nu)interpjac1(nt,ts,nu,da,db,R_or_N);
    [U1,V1] = lowrank(nts,JTM1,ts1,nu,tol,tR,mR);
    %U1 = repmat(sqrt(wghts),1,size(U1,2)).*U1;
    V1 = conj(V1);
    U2=U1;
    U3=U1;
    V2=V1;
    V3=V1;
    %JTM2 = @(ts,nu)interpjac1(nt,ts,nu,da,db,R_or_N);
    %[U2,V2] = lowrank(nts,JTM2,ts2,nu,tol,tR,mR);
    %U2 = diag(sqrt(wghts))*U2;
    %V2 = conj(V2);
    %JTM3 = @(ts,nu)interpjac1(nt,ts,nu,da,db,R_or_N);
    %[U3,V3] = lowrank(nts,JTM1,ts3,nu,tol,tR,mR);
    %U3 = diag(sqrt(wghts))*U3;
    %V3 = conj(V3);
elseif 0 <= opt && opt<1
    [U1,V1] = ID_Cheby(nts,ts1,nu,da,db,tol,1,R_or_N,tR,mR);
    U2=U1;
    U3=U1;
    V2=V1;
    V3=V1;
    %[U2,V2] = ID_Cheby(nts,ts2,nu,da,db,tol,1,R_or_N,tR,mR);
    %[U3,V3] = ID_Cheby(nts,ts3,nu,da,db,tol,1,R_or_N,tR,mR);
    %U = diag(sqrt(wghts))*U;
elseif opt < 0
    [U1,V1] = ID_Cheby(nts,ts1,nu,da,db,tol,-1,R_or_N,tR,mR);
    U2=U1;
    U3=U1;
    V2=V1;
    V3=V1;
    %[U2,V2] = ID_Cheby(nts,ts2,nu,da,db,tol,-1,R_or_N,tR,mR);
    %[U3,V3] = ID_Cheby(nts,ts3,nu,da,db,tol,-1,R_or_N,tR,mR);
    %U = diag(sqrt(wghts))*U;
end
rank1 = size(U1,2);
rank2 = size(U2,2);
rank3 = size(U3,2);
V1 = [zeros(it,rank1);V1];
V2 = [zeros(it,rank2);V2];
V3 = [zeros(it,rank3);V3];
% UU23 = kron(U2,U3);
% VV23 = kron(V2,V3);
% UU12 = kron(U1,U2);
% VV12 = kron(V1,V2);
% UUU = kron(kron(U1,U2),U3);
% VVV = kron(kron(V1,V2),V3);
rank = rank1*rank2*rank3;

wghts1 = wghts;
wghts2 = wghts;
wghts3 = wghts;
vals1 = jacrecur(nts,ts1,wghts1,it-1,da,db).';
% vals1 = vals1*diag(sqrt(wghts));
vals2 = jacrecur(nts,ts2,wghts2,it-1,da,db).';
% vals2 = vals2*diag(sqrt(wghts));
vals3 = jacrecur(nts,ts3,wghts3,it-1,da,db).';
% vals3 = vals3*diag(sqrt(wghts));
% vals12 = kron(vals1,vals2);
% vals23 = kron(vals2,vals3);
% vals = kron(kron(vals1,vals2),vals3);

% Ints = sparse(eye(nts,nts));
% Iit = sparse(eye(it,it));
xs1 = mod(floor(ts1*nts/2/pi),nts)+1;
% Sub1 = sparse(Ints(xs1,:));
xs2 = mod(floor(ts2*nts/2/pi),nts)+1;
% Sub2 = sparse(Ints(xs2,:));
xs3 = mod(floor(ts3*nts/2/pi),nts)+1;

Y = [1:nts]';
X = zeros(nts,1);
S = ones(nts,1);
for j = 1:nts
    X(j) = xs1(j);   
end
P1 = sparse(X,Y,S,nts,nts,nts);
P2 = P1;
P3 = P1;
% PP23 = kron(P2,P3);
% PP12 = kron(P1,P2);
% PPP = kron(PP12,P3);

VV3 = reshape(repmat(V3,nts^2,1),nts,nts^2*rank3);
UU3 = reshape(repmat(U3,nts^2,1),nts,nts^2*rank3);
VV2 = reshape(repmat(V2,nts^2,1),nts,nts^2*rank2);
UU2 = reshape(repmat(U2,nts^2,1),nts,nts^2*rank2);
VV1 = reshape(repmat(V1,nts^2,1),nts,nts^2*rank1);
UU1 = reshape(repmat(U1,nts^2,1),nts,nts^2*rank1);

if  R_or_N > 0
    fun = @(c)invJacPT3d1(c);
else
    ex = exp(1i*nts/2*ts);
    U = U.*repmat(ex,1,rank);
    fun = @(c)invJacPT3d2(c);
end
    
    function y = invJacPT2d(c)
        c = reshape(c,nts,nts^2);
        y = zeros(nts,nts^2);
        y(1:it,:) = vals2*c;

        
        cc = repmat(c,1,rank3);
        dd = UU3.*cc;
        dd = P3*dd;
        fftc = conj(fft(conj(dd)));
        y3 = VV3.*fftc;
        y3 = reshape(y3,nts^3,rank3);
        y3 = reshape(sum(y3,2),nts,nts^2);
        y = y + real(y3);
        for i = 1:nts
            c(:,(i-1)*nts+1:i*nts) = y(:,(i-1)*nts+1:i*nts).';
        end
        
        y = zeros(nts,nts^2);
        y(1:it,:) = vals2*c;

        
        cc = repmat(c,1,rank2);
        dd = UU2.*cc;
        dd = P2*dd;
        fftc = conj(fft(conj(dd)));
        y2 = VV2.*fftc;
        y2 = reshape(y2,nts^3,rank2);
        y2 = reshape(sum(y2,2),nts,nts^2);
        y = y + real(y2);
        for i = 1:nts
            y(:,(i-1)*nts+1:i*nts) = y(:,(i-1)*nts+1:i*nts).';
        end

    end    

    
    function y = invJacPT3d1(c)

	    c = reshape(kron(wghts,kron(wghts,wghts)).*c(:),nts,nts,nts);
        
        c = invJacPT2d(c);
        c = reshape(c,nts^2,nts);
	    c = c.';
        
        y = zeros(nts,nts^2);
        y(1:it,:) = vals1*c;

        cc = repmat(c,1,rank1);
        dd = UU1.*cc;
        dd = P1*dd;
        fftc = conj(fft(conj(dd)));
        y1 = VV1.*fftc;
        y1 = reshape(y1,nts^3,rank1);
        y1 = reshape(sum(y1,2),nts,nts^2);
        y = y + real(y1);
        y = y.';
        y = y(:);
        
    end

end
