function x = BF_idx2vec(siz,iter)
% This code is for indecing.
%
% Copyright 2017 by Yingzhou Li and Haizhao Yang

n = length(siz);
x = zeros(1,n);

%assert(iter <= siz(end),'The iter is out of the range');

for i = 1:n-1
    x(i) = mod(iter-1,siz(i))+1;
    iter = floor((iter-0.5)/siz(i))+1;
end
x(n) = mod(iter-1,siz(n))+1;

end