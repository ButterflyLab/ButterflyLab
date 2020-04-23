% This code computes the inverse matvec for the matrix from the associated
% Legendre polynomials using the butterfly factorization.
%
% Note that the associated Legendre polynomial is even if degree+order is
% even; it's odd if degree+order is odd.

function [a, FacT, AppT] = BF_ALegendre_inv(m,M,g,theta,weight,r,mR,tol,...
    nn,rpt,rand_or_cheb)
% Input:
% m - the order of the matrix slice
% M - 2*N-1 (degree k has a range from m to M)
% g - results of the Legendre transform
% r - the numerical rank of low-rank matrix approximation in the IDBF
% mR - the maximum rank of low-rank approximation by randomized sampling
% tol - set up accuracy
% nn - smalles number of points in leaves
% rpt - # of experiments for application
% rand_or_cheb - whether use random sampling or Mock-Chebshev points for ID
% 
% Output:
% a - the coefficients of the spherical harmonics; a is a cell of length
% 4N-1 for m=-2N+1 to 2N-1; for each m, a{m} is a vector of length 2N-|m|
% for k=|m| to 2N-1
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
ALegendremask = @(k,x,kI,xI) (repmat(xx(xI(x)),1,length(k)) ...
    > real(asin(sqrt(m^2-1/4)./...
    (repmat(kk(kI(k)).',length(x),1)+1/2))))';

G = [g(end/2:-1:1,:), g(1+end/2:end,:)]; gn = size(g, 2);
G1 = [eye(gn); eye(gn)]; G2 = [-eye(gn);eye(gn)];
Lodd = length(kOddIdx); Leven = length(kEvenIdx);
N = numel(theta) / 2;
xIdx = (N+1) : (N*2);

% odd
fun1 = @(k,x) fun(xIdx(x),kOddIdx(k))';
mask = @(k,x) ALegendremask(k,x,kOddIdx,xIdx);
t = tic;
FOdd = IDBF_mask(fun1,kk(kOddIdx),xx(xIdx),mask,nn,r,mR,tol,rand_or_cheb);
FacT = FacT + toc(t);
for i = 1 : rpt
    t = tic;
    aOdd = BF_mask_apply(FOdd, G * G1, Lodd);
    AppT(i) = AppT(i) + toc(t);
end

% even
if M>m
    fun1 = @(k,x) fun(xIdx(x),kEvenIdx(k))';
    mask = @(k,x) ALegendremask(k,x,kEvenIdx,xIdx);
    t = tic;
    FEven = IDBF_mask(fun1,kk(kEvenIdx),xx(xIdx),mask,nn,r,mR,tol,rand_or_cheb);
    FacT = FacT + toc(t);
    for i = 1 : rpt
        t = tic;
        aEven = BF_mask_apply(FEven, G * G2, Leven);
        AppT(i) = AppT(i) + toc(t);
    end
end

% combine
a = zeros(M-m+1,size(g,2));
a(1:2:end,:) = aOdd;
if M>m
    a(2:2:end,:) = aEven;
end

end
