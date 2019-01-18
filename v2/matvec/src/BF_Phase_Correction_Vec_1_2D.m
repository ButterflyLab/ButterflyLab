function [U,disPosu] = BF_Phase_Correction_Vec_1_2D(U,disPosu,tau,thre)
% O(N log N) operation and memory complexity.
% Reference: H. Yang, A Unified Framework for Oscillatory Integral
% Transform: When to use NUFFT or Butterfly factorization? preprint, 2018.

if nargin < 4, thre = 1e10; end
pos = find(disPosu(:,2)>=2);
numB = numel(pos);
for cnt = 1:numB
    st = pos(cnt);
    if cnt < numB
        ed = pos(cnt+1)-1;
    else
        ed = length(U);
    end
    tmp = U(st:ed);
    len = disPosu(st+1,1)-disPosu(st,1);
    tmp = reshape(tmp,len);
    [posu,posv] = size(tmp);
    posu = 1:posu; posv = 1:posv;
    tmp = BF_Phase_Correction_Mat_1D(tmp,tmp.',posu,posv,1,1,1,1,tau);
    U(st:ed) = tmp(:);
end
end







