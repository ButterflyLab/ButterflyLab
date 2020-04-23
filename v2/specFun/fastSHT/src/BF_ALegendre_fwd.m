% This code computes the matvec for the matrix from the associated Legendre
% polynomials using the butterfly factorization.
%
% Note that the associated Legendre polynomial is even if degree+order is
% even; it's odd if degree+order is odd.

function [g, FacT, AppT] = BF_ALegendre_fwd(m,M,a,theta,weight,r,mR,tol,...
    nn,rpt,rand_or_cheb)
% Input:
% m - the order of the matrix slice
% M - 2*N-1 (degree k has a range from m to M)
% a - the coefficients of the spherical harmonics; a is a cell of length
% 4N-1 for m=-2N+1 to 2N-1; for each m, a{m} is a vector of length 2N-|m|
% for k=|m| to 2N-1
% r - the numerical rank of low-rank matrix approximation in the IDBF
% mR - the maximum rank of low-rank approximation by randomized sampling
% tol - set up accuracy
% nn - smalles number of points in leaves
% rpt - # of experiments for application
% rand_or_cheb - whether use random sampling or Mock-Chebshev points for ID
% 
% Output:
% g - results of the Legendre transform
% FacT - factorization time
% AppT - application time

if nargin < 8, tol = 1e-10; end
if nargin < 9, nn = 256; end
if nargin < 10, rpt = 1; end
if nargin < 11, rand_or_cheb = 'cheb'; end

FacT = 0;
AppT = zeros(1, rpt);

kk = (m:M)'; kOddIdx = 1:2:(M-m+1); kEvenIdx = 2:2:(M-m+1);
xx = theta;
fun = @(x,k) ALegendrefun1(x,k,m,xx,kk,weight);
ALegendremask = @(x,k,xI,kI) repmat(xx(xI(x)),1,length(k)) ...
    > real(asin(sqrt(m^2-1/4)./...
    (repmat(kk(kI(k)).',length(x),1)+1/2)));

N = numel(theta) / 2;
xIdx = (N+1) : (N*2);

% odd
fun1 = @(x,k) fun(xIdx(x),kOddIdx(k));
mask = @(x,k) ALegendremask(x,k,xIdx,kOddIdx);
t = tic;
FOdd = IDBF_mask(fun1,xx(xIdx),kk(kOddIdx),mask,nn,r,mR,tol,rand_or_cheb);
FacT = FacT + toc(t);
for i = 1 : rpt
    t = tic;
    gOdd = BF_mask_apply(FOdd, a(1:2:end,:), N);
    AppT(i) = AppT(i) + toc(t);
end

% even
if M>m
    fun1 = @(x,k) fun(xIdx(x),kEvenIdx(k));
    mask = @(x,k) ALegendremask(x,k,xIdx,kEvenIdx);
    t = tic;
    FEven = IDBF_mask(fun1,xx(xIdx),kk(kEvenIdx),mask,nn,r,mR,tol,rand_or_cheb);
    FacT = FacT + toc(t);
    for i = 1 : rpt
        t = tic;
        gEven = BF_mask_apply(FEven, a(2:2:end,:), N);
        AppT(i) = AppT(i) + toc(t);
    end
end

% combine
g = zeros(2*N,size(a,2));
if M>m
    g(end/2+1:end,:) = gOdd+gEven;
    g(1:end/2,:) = gOdd(end:-1:1,:)-gEven(end:-1:1,:);
else
    g(end/2+1:end,:) = gOdd;
    g(1:end/2,:) = gOdd(end:-1:1,:);
end

end
