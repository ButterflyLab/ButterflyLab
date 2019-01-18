function cidxs = BF_childidx(siz,idx)
% This code computes the children indecies
%
% Copyright 2017 by Yingzhou Li and Haizhao Yang

n = length(siz);
cidxs = zeros(1,2^n);
for itc = 1:2^n
    vec_child = BF_idx2vec(2*ones(1,n),itc);
    vec = BF_idx2vec(siz,idx);
    cidxs(itc) = BF_vec2idx(2*siz, (vec-1)*2 + vec_child);
end

end