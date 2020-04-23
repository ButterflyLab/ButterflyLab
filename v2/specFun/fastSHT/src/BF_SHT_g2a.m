function [a, FacT, AppT] = BF_SHT_g2a(g,N,r,mR,tol,nn,rpt)
% computing a(k,m) from g(theta,m) using butterfly factorization
% g(theta,m)=\sum_{k=|m|}^{2N-1}a(k,m)\tilte{P}_k^{|m|}(cos(\theta))
%
% Input:
% g - results of the Legendre transform
% N - the degree of the spherical harmonic transform is truncated at -2N+1
% and 2N-1
% r - the numerical rank of low-rank matrix approximation in the IDBF
% mR - the maximum rank of low-rank approximation by randomized sampling
% tol - set up accuracy
% nn - smalles number of points in leaves
% rpt - # of experiments for application
% 
% Output:
% a - the coefficients of the spherical harmonics; a is a cell of length
% 4N-1 for m=-2N+1 to 2N-1; for each m, a{m} is a vector of length 2N-|m|
% for k=|m| to 2N-1
% FacT - factorization time
% AppT - application time

if nargin < 5, tol = 1e-10; end
if nargin < 6, nn = 256; end
if nargin < 7, rpt = 1; end

[theta,weight] = legenQuad(2*N); 
theta = acos(real(theta(:))); weight = real(weight(:));

a = cell(4*N-1,1);
for m = -2*N+1:2*N-1
    a{m+2*N} = zeros(2*N-abs(m),1);
end

A = cell(2*N, 1); ft = zeros(2*N, 1); at = zeros(2*N, rpt);
for m = 1:2*N-1 % m is the order
    [A{m}, ft(m), at(m,:)] = ...
        BF_ALegendre_inv(abs(m),2*N-1,[g(:,m+2*N),g(:,-m+2*N)],...
        theta,weight,r,mR,tol,nn,rpt);
end
[A{2*N}, ft(2*N), at(2*N,:)] = ...
    BF_ALegendre_inv(0,2*N-1,g(:,2*N),theta,weight,r,mR,tol,nn,rpt);

for m = 1:2*N-1
    a{m+2*N} = A{m}(:,1); a{-m+2*N} = A{m}(:,2);
end
a{2*N} = A{2*N};

FacT = sum(ft); AppT = sum(at);

end