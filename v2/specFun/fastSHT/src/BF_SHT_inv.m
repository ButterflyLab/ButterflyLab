% This code computes the inverse spherical harmonic transform using the
% butterfly factorization.
%
% The spherical harmonic transform is defined by
% f(psi,theta) = \sum_{k=0}^{2N-1} \sum_{m=-k}^{k} a(k,m)
% \tilte{P}_k^{|m|}(cos(theta)) exp(i*m*psi)
%
% Forward transform computes a(k,m) from f(psi,theta)
% Inverse transform computes f(psi,theta) from a(k,m)
%
% The spherical harmonic transform can be written as
% f(psi,theta) = \sum_{m=-2N+1}^{2N-1} g(theta,m)exp(i*m*psi),
% where
% g(theta,m)=\sum_{k=|m|}^{2N-1}a(k,m)\tilte{P}_k^{|m|}(cos(\theta))
%
% theta_i, i=0,...,2*N-1 are the 2N roots of \tilte{P}_{2N}^0(cos(theta_i))=0.
% psi_i, i=0,...,4*N-2  = 2*pi*(i+1/2)/(4*N-1)

function f = BF_SHT_inv(a,N,r,mR,tol)
%
% Input:
% a - the coefficients of the spherical harmonics; a is a cell of length
% 4N-1 for m=-2N+1 to 2N-1; for each m, a{m} is a vector of length 2N-|m|
% for k=|m| to 2N-1
% N - the degree of the spherical harmonic transform is truncated at -2N+1
% and 2N-1
% r - the numerical rank of low-rank matrix approximation in the IDBF
% mR - the maximum rank for low-rank approximation by randomized sampling
% tol - set up accuracy
%
% Output:
% f - the reconstructed function f(psi,theta), a matrix of size 4*N-1 by 2*N

assert(length(a)==4*N-1,'the size of coefficients does not match');

% reconstructing g(theta,m) from a(k,m) using butterfly factorization
% g(theta,m)=\sum_{k=|m|}^{2N-1}a(k,m)\tilte{P}_k^{|m|}(cos(\theta))
g = BF_SHT_a2g(a,N,r,mR,tol);

% reconstructing f(psi,theta) from g(theta,m) using FFT
% f(psi,theta) = \sum_{m=-2N+1}^{2N-1} g(theta,m)exp(i*m*psi)
f = BF_SHT_g2f(g,N);

end
