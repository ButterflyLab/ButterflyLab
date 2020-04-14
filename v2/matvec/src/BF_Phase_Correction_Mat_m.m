% Multidimensional Phase Recovery and Interpolative Decomposition 
% Butterfly Factorization, preprint, 2019.
% Copyright 2018 by Ze Chen, Juan Zhang, Ken Ho, and Haizhao Yang

% O(N log N) operation and memory complexity.

function [M1,M2] = BF_Phase_Correction_Mat_m(M1,M2,cs,rs,P1,P1c,P2,P2c,tau)
% M2 for rows and M1 for columns, all talk skinny matrices

    % first row
    M2(:,1) = BF_Phase_Correction_Vec_m(M2(:,1),P2,tau);
    % first column
    M1(:,1) = BF_Phase_Correction_Vec_m(M1(:,1),P1,tau);
    % the rest rows
    M2(1, 2:end) = M1(rs(2:end), 1)';
    for ri = 2 : length(rs)
        M2(:,ri) = BF_Phase_Correction_Vec_m(M2(:,ri),P2,tau);
    end
    % the rest columns
    M1(rs, 2:end) = M2(cs(2:end), :)';
    for ci = 2 : length(cs)
        M1(:,ci) = BF_Phase_Correction_Vec_m(M1(:,ci),P1c,tau);
    end
    
end