% This code computes the forward spherical harmonic transform using the
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

function a = BF_SHT_fwd(f,N,r,mR,tol)
%
% Input:
% f - the reconstructed function f(psi,theta), a matrix of size 4*N-1 by 2*N
% N - the degree of the spherical harmonic transform is truncated at -2N+1
% and 2N-1
% r - the numerical rank of low-rank matrix approximation in the IDBF
% mR - the maximum rank for low-rank approximation by randomized sampling
% tol - set up accuracy
%
% Output:
% a - the coefficients of the spherical harmonics; a is a cell of length
% 4N-1 for m=-2N+1 to 2N-1; for each m, a{m} is a vector of length 2N-|m|
% for k=|m| to 2N-1

assert(norm(size(f)-[4*N-1,2*N])<1e-5,'the size of the given function does not match');

% compute g(theta,m) from f(psi,theta) using FFT
% f(psi,theta) = \sum_{m=-2N+1}^{2N-1} g(theta,m)exp(i*m*psi)
g = BF_SHT_f2g(f,N);

% computing a(k,m) from g(theta,m) using butterfly factorization
% g(theta,m)=\sum_{k=|m|}^{2N-1}a(k,m)\tilte{P}_k^{|m|}(cos(\theta))
a = BF_SHT_g2a(g,N,r,mR,tol);

end
