% Multidimensional Phase Recovery and Interpolative Decomposition 
% Butterfly Factorization, preprint, 2019.
% Copyright 2018 by Ze Chen, Juan Zhang, Ken Ho, and Haizhao Yang

% O(N log N) operation and memory complexity.

function U = BF_Phase_Correction_Vec_m(U,P,disPosu,tau)
% U - phase vector for recovery
% P - recovery path matrix
% disPosu - discontinuity locations of U
% tau - scale value of phase vector

N = numel(U);

for c = 1 : N-1
    bg = P(c,1); ed = P(c,2);
    if ~ismember(ed,disPosu)
        U(ed) = U(ed) - round((U(ed)-U(bg))/tau)*tau;
    end
end

end