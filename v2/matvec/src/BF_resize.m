function Factor = BF_resize(Factor, N1, N2)
% This code resizes the dimension of the factors in BF to make sure consistency.
% the BF matrix is of size N1 by N2.
%
% Copyright 2018 Haizhao Yang

if nargin< 3, N2 = N1; end

sz = size(Factor.V,2);
if sz < N2
    Factor.V(1,N2) = 0;
end
N2 = size(Factor.V,1);

for i=length(Factor.HTol):-1:1
    sz = size(Factor.HTol{i},2);
    if sz < N2
        Factor.HTol{i}(1,N2) = 0;
    else if sz > N2
            if i < length(Factor.HTol)
                Factor.HTol{i+1}(sz,1) = 0;
            else
                Factor.V(sz,1) = 0;
            end
        end
    end
    N2 = size(Factor.HTol{i},1);
end

sz = size(Factor.M,2);
if sz < N2
    Factor.M(1,N2) = 0;
else if sz > N2
        Factor.HTol{1}(sz,1) = 0;
    end
end
N2 = size(Factor.M,1);

for i=1:length(Factor.GTol)
    sz = size(Factor.GTol{i},2);
    if sz < N2
        Factor.GTol{i}(1,N2) = 0;
    else if sz > N2
            if i > 1
                Factor.GTol{i-1}(sz,1) = 0;
            else
                Factor.M(sz,1) = 0;
            end
        end
    end
    N2 = size(Factor.GTol{i},1);
end

sz = size(Factor.U,2);
if sz < N2
    Factor.U(1,N2) = 0;
else if sz > N2
        if length(Factor.GTol)>0
            Factor.GTol{end}(sz,1) = 0;
        else
            Factor.M(sz,1) = 0;
        end
    end
end
if size(Factor.U,1) < N1
    Factor.U(N1,1) = 0;
end

end
