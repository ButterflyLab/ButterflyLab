function  y = HSSBF_mat_apply(Factor,p)
% HSSBF_apply applies a HSSBF matrix stored in a sparse matrix 
% factorization format in Factor to a give vector p.
%
% Input:
% p: a given vector
% Factor: a structure storing the matrix factors, i.e.
%       A \approx Factor{1}*...*Factor{n}
%
% Output:
% y: A*p
%
%
% Reference: A Hierarchical Butterfly LU Preconditioner for Two-Dimensional
% Electromagnetic Scattering Problems Involving Open Surfaces, preprint,
% 2019
% Copyright reserved.
% By Haizhao Yang, 2018, National University of Singapore

n = length(Factor);
y = p;
for cnt = n:-1:1
    y = Factor{cnt}*y;
end