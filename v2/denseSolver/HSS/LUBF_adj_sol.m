function y = LUBF_adj_sol(Factor, p,L_or_U)
% LUBF_adj_sol applies the data-sparse inverse of the conjugate transpose
% of a triangular HSSBF matrix to a vector p
%
% Input:
% Factor: a structure storing the data-sparse representation of a lower or
%         upper triangular matrix in a HSSBF format, and some auxiliary
%         variables.
% Factor.isLeaf: 1, if all matrices are dense
%                0, if all matrices are data-sparse
% Factor.sz: the size of the matrix stored in Factor
% Factor.szsub: the size of the matrix stored in A11
% Factor.A11
% Factor.A22
% Factor.A12
% Factor.A21
% p: right hand side vector.
% L_or_U: indicate the matrix in Factor is a lower or upper triangular
% matrix
%
% Output:
% y - inv(Factor')*p
%
%
% Reference: A Hierarchical Butterfly LU Preconditioner for Two-Dimensional
% Electromagnetic Scattering Problems Involving Open Surfaces, preprint,
% 2019
% Copyright reserved.
% By Haizhao Yang, 2018, National University of Singapore

n1 = Factor.szsub;
switch L_or_U
    case 'U'
        if Factor.isLeaf
            x1 = (Factor.A11')\p(1:n1,:);
            x2 = (Factor.A22')\(p((n1+1):end,:)-Factor.A12'*x1);
            y = [x1;x2];
        else
            x1 = LUBF_adj_sol(Factor.A11,p(1:n1,:),L_or_U);
            x2 = LUBF_adj_sol(Factor.A22,p((n1+1):end,:)-BF_adj_apply(Factor.A12,x1),L_or_U);
            y = [x1;x2];
        end
    case 'L'
        if Factor.isLeaf
            x2 = Factor.A22'\p((n1+1):end,:);
            x1 = Factor.A11'\(p(1:n1,:)-Factor.A21'*x2);
            y = [x1;x2];
        else
            x2 = LUBF_adj_sol(Factor.A22,p((n1+1):end,:),L_or_U);
            x1 = LUBF_adj_sol(Factor.A11,p(1:n1,:)-BF_adj_apply(Factor.A21,x2),L_or_U);
            y = [x1;x2];
        end
end
end
