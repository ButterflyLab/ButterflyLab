function idx = vec2idx(siz,vec)
% This code is for indecing.
%
% Copyright 2017 by Yingzhou Li and Haizhao Yang

%assert(all(vec <= siz),'The iter is out of the range');

n = length(siz);
siz = cumprod(siz);

idx = vec(1);
for i = 2:n
    idx = idx + (vec(i)-1)*siz(i-1);
end

end