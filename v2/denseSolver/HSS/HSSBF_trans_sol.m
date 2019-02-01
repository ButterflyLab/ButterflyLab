function y = HSSBF_trans_sol(Factor, p)
% HSSBF_trans_sol applies the transpose of the data-sparse inverse of a 
% HSSBF matrix to a vector p.
%
% Input:
% Factor: a structure storing the data-sparse representation of the inverse
%         of a HSSBF matrix A, and some auxiliary variables.
% Factor.isLeaf: 1, if all matrices are dense
%                0, if all matrices are data-sparse
% Factor.sz: the size of the matrix stored in Factor
% Factor.B
% Factor.Btr
% Factor.Bctr
% Factor.A11
% Factor.A22
% Factor.A11inv
% Factor.A22inv
% p: right hand side vector.
%
% Output:
% y - inv(A).'*p
%
%
% Reference: A Hierarchical Butterfly LU Preconditioner for Two-Dimensional
% Electromagnetic Scattering Problems Involving Open Surfaces, preprint,
% 2019
% Copyright reserved.
% By Haizhao Yang, 2018, National University of Singapore

n1 = Factor.A11inv.sz;
if Factor.isLeaf
    p = (Factor.B.')*p;
    A11invp = (Factor.A11.')\p(1:n1,:);
    A22invp = (Factor.A22.')\p((n1+1):end,:);
    y = [A11invp;A22invp];
else
    p = Factor.Btr(p);
    A11invp = HSSBF_trans_sol(Factor.A11inv,p(1:n1,:));
    A22invp = HSSBF_trans_sol(Factor.A22inv,p((n1+1):end,:));
    y = [A11invp;A22invp];
end
end
