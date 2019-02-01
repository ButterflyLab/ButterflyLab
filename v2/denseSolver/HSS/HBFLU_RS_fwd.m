function  [L,U] = HBFLU_RS_fwd(Afun,x,k,rk,tol,lsz,method)
% HBFLU_RS_fwd compresses two L and U such that L+U-I=A, where I is an
% identity, L is a lower triangular matrix, U is a upper triangular matrix, 
% and A is given by Afun.
%
% Input:
% Afun: a function handle for evaluating an arbitrary entry of the HBFLU
%       matrix A in O(1) operations, i.e. A(i,j) = Afun(i,j).
% x and k: denote the row and column indices of the submatrix of A to be
%       inverted and compressed, i.e. invert and compress the submatrix
%       A(x,k).
% rk:    the maximum rank of the low-rank approximation
% tol:  the accuracy parameter of the low-rank approximation. If we can
%       achieve a low-rank approximation with a rank smaller than rk and an
%       accuracy tol, we will use a smaller rank.
% lsz:  a size parameter; when the matrix A has a size <= lsz, we
%       don't compress
% method: specify which BF to computer HBFLU
%         1: IDBF
%         2: CURBF
%         3: CURBF2
%
% Output:
% L: a structure storing the data-sparse representation of the matrix
%        tril(A), and some auxiliary variables.
% L.isLeaf: 1, if all matrices are dense
%                0, if all matrices are data-sparse
% L.sz: the size of A
% L.szsub: the size of A11
% L.A11
% L.A12
% L.A21
% L.A22
% U: a structure storing the data-sparse representation of the matrix
%        triu(A)-diag(A)+I, and some auxiliary variables.
% U.isLeaf: 1, if all matrices are dense
%                0, if all matrices are data-sparse
% U.sz: the size of A
% U.szsub: the size of A11
% U.A11
% U.A12
% U.A21
% U.A22
%
%
% Reference: A Hierarchical Butterfly LU Preconditioner for Two-Dimensional
% Electromagnetic Scattering Problems Involving Open Surfaces, preprint,
% 2019
% Copyright reserved.
% By Haizhao Yang, 2018, National University of Singapore

if nargin < 7, method = 1; end

M = numel(x); N = numel(k);
L = struct('isLeaf',[],'sz',[],'szsub',[],'A11',[],'A12',[],'A21',[],'A22',[]);
U = struct('isLeaf',[],'sz',[],'szsub',[],'A11',[],'A12',[],'A21',[],'A22',[]);
if M~=N
    error('matrix is not square');
end
if M>lsz
    L.sz = M;
    L.szsub = numel(x(1:end/2));
    L.isLeaf = 0;
    U.sz = N;
    U.szsub = numel(x(1:end/2));
    U.isLeaf = 0;
    [L.A11,U.A11] = HBFLU_RS_fwd(Afun,x(1:end/2),k(1:end/2),rk,tol,lsz,method);
    [L.A22,U.A22] = HBFLU_RS_fwd(Afun,x((end/2+1):end),k((end/2+1):end),rk,tol,lsz,method);
    switch method
        case 1
            U.A12 = BF_IDBF(Afun,x(1:end/2),k((end/2+1):end),8,rk,tol);
            L.A21 = BF_IDBF(Afun,x((end/2+1):end),k(1:end/2),8,rk,tol);
        case 2
            U.A12 = CURBF(Afun,x(1:end/2),k((end/2+1):end),rk,tol,0);
            L.A21 = CURBF(Afun,x((end/2+1):end),k(1:end/2),rk,tol,0);
        case 3
            U.A12 = CURBF2(Afun,x(1:end/2),k((end/2+1):end),max(rk,8),rk,tol);
            L.A21 = CURBF2(Afun,x((end/2+1):end),k(1:end/2),max(rk,8),rk,tol);
    end
else
    L.sz = M;
    L.szsub = numel(x(1:end/2));
    L.isLeaf = 1;
    U.sz = N;
    U.szsub = numel(x(1:end/2));
    U.isLeaf = 1;
    
    temp = Afun(x(1:end/2),k(1:end/2));
    L.A11 = tril(temp);
    U.A11 = triu(temp)-diag(diag(temp)-ones(numel(x(1:end/2)),1));
    
    L.A21 = Afun(x((end/2+1):end),k(1:end/2));
    U.A12 = Afun(x(1:end/2),k((end/2+1):end));
    
    temp = Afun(x((end/2+1):end),k((end/2+1):end));
    L.A22 = tril(temp);
    U.A22 = triu(temp)-diag(diag(temp)-ones(numel(x((end/2+1):end)),1));
end
