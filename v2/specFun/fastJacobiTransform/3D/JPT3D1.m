function [fun,rank1,rank2,rank3] = JPT3D1(nts,da,db,tR,mR,tol,opt,R_or_N)
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
%  fun    -  a function handle computing the 3D forward uniform Jacobi polynomial transform such that
%                 fun(c(:)) = kron(kron(J,J),J)*c(:), 
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
    %U = diag(sqrt(wghts))*U;
    V1 = conj(V1);
    U2=U1;
    U3=U1;
    V2=V1;
    V3=V1;
    %JTM2 = @(ts,nu)interpjac1(nt,ts,nu,da,db,R_or_N);
    %[U2,V2] = lowrank(nts,JTM2,ts2,nu,tol,tR,mR);
    %U = diag(sqrt(wghts))*U;
    %V2 = conj(V2);
    %JTM3 = @(ts,nu)interpjac1(nt,ts,nu,da,db,R_or_N);
    %[U3,V3] = lowrank(nts,JTM1,ts3,nu,tol,tR,mR);
    %U = diag(sqrt(wghts))*U;
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
V3 = [zeros(it,rank2);V3];
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
vals1 = jacrecur(nts,ts1,wghts1,it-1,da,db);
vals2 = jacrecur(nts,ts2,wghts2,it-1,da,db);
vals3 = jacrecur(nts,ts3,wghts3,it-1,da,db);
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
% Sub3 = sparse(Ints(xs3,:));
% SSS = kron(kron(Sub1,Sub2),Sub3);


% xsub23 = zeros(nts*nts,1);
% for i = 1:nts
%     for j =1:nts
%         xsub23((i-1)*nts+j) = (xs2(i)-1)*nts + xs3(j);
%     end
% end
% xsub12 = zeros(nts*nts,1);
% for i = 1:nts
%     for j =1:nts
%         xsub12((i-1)*nts+j) = (xs1(i)-1)*nts + xs2(j);
%     end
% end

VV3 = reshape(repmat(V3,nts^2,1),nts,nts^2*rank3);
UU3 = reshape(repmat(U3,nts^2,1),nts,nts^2*rank3);
VV2 = reshape(repmat(V2,nts^2,1),nts,nts^2*rank2);
UU2 = reshape(repmat(U2,nts^2,1),nts,nts^2*rank2);
VV1 = reshape(repmat(V1,nts^2,1),nts,nts^2*rank1);
UU1 = reshape(repmat(U1,nts^2,1),nts,nts^2*rank1);



if  R_or_N > 0
    fun = @(c)JacPT3d1(c);
else
    %ex = exp(1i*nts/2*ts);
    %U = U.*repmat(ex,1,rank);
    fun = @(c)JacPT3d1(c);
end
   
   function y = JacPT2d(c)
        c = reshape(c,nts,nts^2);
        y = vals3*c(1:it,:);

        cc = repmat(c,1,rank3);
        dd = VV3.*cc;
        fftc = nts*ifft(dd);
        fftc = fftc(xs3,:);
        y3 = UU3.*fftc;
        y3 = reshape(y3,nts^3,rank3);
        y3 = reshape(sum(y3,2),nts,nts^2);
        y = y + real(y3);
        for i = 1:nts
            c(:,(i-1)*nts+1:i*nts) = y(:,(i-1)*nts+1:i*nts).';
        end
        
        y = vals2*c(1:it,:);
        
        cc = repmat(c,1,rank2);
        dd = VV2.*cc;

        fftc = nts*ifft(dd);
        fftc = fftc(xs2,:);
        y2 = UU2.*fftc;
        y2 = reshape(y2,nts^3,rank2);
        y2 = reshape(sum(y2,2),nts,nts^2);
        y = y + real(y2);
        for i = 1:nts
            y(:,(i-1)*nts+1:i*nts) = y(:,(i-1)*nts+1:i*nts).';
        end

    end

    function y = JacPT3d1(c)


        c = JacPT2d(c);
        c = reshape(c,nts^2,nts);
	    c = c.';

        y = vals1*c(1:it,:);
        cc = repmat(c,1,rank1);
        dd = VV1.*cc;
        fftc = nts*ifft(dd);
        fftc = fftc(xs1,:);
        y1 = UU1.*fftc;
        y1 = reshape(y1,nts^3,rank1);
        y1 = reshape(sum(y1,2),nts,nts^2);
        y = y + real(y1);
        y = y.';
        y = y(:);
    end
    

    
end
