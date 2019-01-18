function pidx = BF_parentidx(siz,idx)
% This code identifies the parent indices.
%
% Copyright 2017 by Yingzhou Li and Haizhao Yang

vec = BF_idx2vec(siz,idx);
pidx = BF_vec2idx(siz/2, floor((vec+1)/2));

end