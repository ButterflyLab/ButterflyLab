function  Factor = HSSBF_MV_fwd(Afun,x,k,r,tol,lsz)
% HSSBF_MV_fwd compresses a HSSBF matrix stored in a function handle Afun and 
% returns a data-sparse HSSBF format in Factor.
%
% Input:
% Afun: a function handle for evaluating an arbitrary entry of the HSSBF
%       matrix A in O(1) operations, i.e. A(i,j) = Afun(i,j).
% x and k: denote the row and column indices of the submatrix of A to be
%       inverted and compressed, i.e. invert and compress the submatrix
%       A(x,k).
% r:    the maximum rank of the low-rank approximation
% tol:  the accuracy parameter of the low-rank approximation. If we can
%       achieve a low-rank approximation with a rank smaller than r and an
%       accuracy tol, we will use a smaller rank.
% lsz:  a size parameter; when the matrix A has a size <= lsz, we
%       don't compress
%
% Output:
% Factor: a structure storing the data-sparse representation of the matrix
%        A(x,k), and some auxiliary variables.
% Factor.isLeaf: 1, if all matrices are dense
%                0, if all matrices are data-sparse
% Factor.sz: the size of A
% Factor.szsub: the size of A11
% Factor.A11
% Factor.A12
% Factor.A21
% Factor.A22
%
%
% Reference: A Hierarchical Butterfly LU Preconditioner for Two-Dimensional
% Electromagnetic Scattering Problems Involving Open Surfaces, preprint,
% 2019
% Copyright reserved.
% By Haizhao Yang, 2018, National University of Singapore

M = numel(x); N = numel(k);
Factor = struct('isLeaf',[],'sz',[],'szsub',[],'A11',[],'A12',[],'A21',[],'A22',[]);
if M~=N
    error('matrix is not square');
end
% The idea doesn't work: how to devide a function handle?
% if M>=lsz
%     Factor.sz = M;
%     Factor.szsub = numel(x(1:end/2));
%     Factor.isLeaf = 0;
%     Factor.A11 = HSSBF_MV_fwd(Afun,x(1:end/2),k(1:end/2),r,tol,lsz);
%     Factor.A22 = HSSBF_MV_fwd(Afun,x((end/2+1):end),k((end/2+1):end),r,tol,lsz);
%     Factor.A12 = CURBF(Afun,x(1:end/2),k((end/2+1):end),r,tol);
%     Factor.A21 = CURBF(Afun,x((end/2+1):end),k(1:end/2),r,tol);
% else
%     Factor.sz = M;
%     Factor.szsub = numel(x(1:end/2));
%     Factor.isLeaf = 1;
%     Factor.A11 = Afun(x(1:end/2),k(1:end/2));
%     Factor.A22 = Afun(x((end/2+1):end),k((end/2+1):end));
%     Factor.A12 = Afun(x(1:end/2),k((end/2+1):end));
%     Factor.A21 = Afun(x((end/2+1):end),k(1:end/2));
% end