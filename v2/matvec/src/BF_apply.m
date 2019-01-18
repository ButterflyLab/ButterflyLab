function y = BF_apply(Factor, p)
% This code apply the butterfly factorization (stored in Factor) to a
% vector p.
%
% Copyright
% 2018 by Qiyuan Pang and Haizhao Yang
% 2017 by Yingzhou Li and Haizhao Yang

switch numel(fieldnames(Factor))
    case 5
        y = Factor.V*p;
        
        for i=length(Factor.HTol):-1:1
            y = Factor.HTol{i}*y;
        end
        
        y = Factor.M*y;
        
        for i=1:length(Factor.GTol)
            y = Factor.GTol{i}*y;
        end
        
        y = Factor.U*y;
        
    case 4
        L = length(Factor.U);
        assert(length(Factor.V) == L)  % check size
        y = p;
        for l = 1:L
            y = Factor.V{l}*y;  % apply V factors
        end
        y = Factor.S*y;         % apply S factor
        for l = L:-1:1
            y = Factor.U{l}*y;  % apply U factors
        end
        
    case 3
        y = p;
        r = size(Factor.U,2);
        for ii = 1:r
            y = Factor.V{ii}*y;
        end
        y = Factor.S*y;
        for ii = r:-1:1
            y = Factor.U{ii}*y;
        end
end

