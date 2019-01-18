function y = BF_adj_apply(Factor, p)
% This code apply the adjoint of the butterfly factorization (stored in
% Factor) to a vector p.
%
% Copyright 2018 by Haizhao Yang
% Copyright 2017 by Yingzhou Li and Haizhao Yang
switch numel(fieldnames(Factor))
    case 5
        y = Factor.U'*p;
        
        for i=length(Factor.GTol):-1:1
            y = Factor.GTol{i}'*y;
        end
        
        y = Factor.M'*y;
        
        for i=1:length(Factor.HTol)
            y = Factor.HTol{i}'*y;
        end
        
        y = Factor.V'*y;
    case 4
        L = length(Factor.U);
        assert(length(Factor.V) == L)  % check size
        y = p;
        for l = 1:L
            y = Factor.U{l}'*y;  % apply U factors
        end
        y = Factor.S'*y;         % apply S factor
        for l = L:-1:1
            y = Factor.V{l}'*y;  % apply V factors
        end
        
    case 3
        y = p;
        r = size(Factor.U,2);
        for ii = 1:r
            y = Factor.U{ii}'*y;
        end
        y = Factor.S'*y;
        for ii = r:-1:1
            y = Factor.V{ii}'*y;
        end
end
end
