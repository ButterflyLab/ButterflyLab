function y = BF_adj_apply(Factor, p)
% This code apply the transpose of the butterfly factorization (stored in 
% Factor) to a vector p.
%
% Copyright 2017 by Yingzhou Li and Haizhao Yang

y = Factor.U.'*p;

for i=length(Factor.GTol):-1:1
    y = Factor.GTol{i}.'*y;
end

y = Factor.M.'*y;

for i=1:length(Factor.HTol)
    y = Factor.HTol{i}.'*y;
end

y = Factor.V.'*y;

end
