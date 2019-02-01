function  y = HSSBF_adj_apply(Factor,p)
% HSSBF_apply applies the adjoint of HSSBF matrix stored in a data-sparse
% HSSBF format in Factor to a give vector p.
%
% Input:
% p: a given vector
% Factor: a structure storing the data-sparse representation of the matrix
%         A(x,k), and some auxiliary variables.
% Factor.isLeaf: 1, if all matrices are dense
%                0, if all matrices are data-sparse
% Factor.sz: the size of A
% Factor.szsub: the size of A11
% Factor.A11
% Factor.A12
% Factor.A21
% Factor.A22
%
% Output:
% y: A'*p
%
%
% Reference: A Hierarchical Butterfly LU Preconditioner for Two-Dimensional
% Electromagnetic Scattering Problems Involving Open Surfaces, preprint,
% 2019
% Copyright reserved.
% By Haizhao Yang, 2018, National University of Singapore

y = zeros(size(p));
M = Factor.sz;
if M~=size(p,1)
    error('matvec dimension problem');
end
M1 = Factor.szsub;
if Factor.isLeaf
    if isempty(Factor.A21) == 1
        y(1:M1,:) = Factor.A11'*p(1:M1,:);
    else
        y(1:M1,:) = Factor.A11'*p(1:M1,:)+Factor.A21'*p((1+M1):end,:);
    end
    if isempty(Factor.A12) == 1
        y((1+M1):end,:) = Factor.A22'*p((1+M1):end,:);
    else
        y((1+M1):end,:) = Factor.A12'*p(1:M1,:)+Factor.A22'*p((1+M1):end,:);
    end
else
    if isempty(Factor.A21) == 1
        y(1:M1,:) = HSSBF_adj_apply(Factor.A11,p(1:M1,:));
    else
        y(1:M1,:) = HSSBF_adj_apply(Factor.A11,p(1:M1,:)) + BF_adj_apply(Factor.A21,p((M1+1):end,:));
    end
    if isempty(Factor.A12) == 1
        y((1+M1):end,:) = HSSBF_adj_apply(Factor.A22,p((1+M1):end,:));
    else
        y((1+M1):end,:) = HSSBF_adj_apply(Factor.A22,p((1+M1):end,:)) + BF_adj_apply(Factor.A12,p(1:M1,:));
    end
end