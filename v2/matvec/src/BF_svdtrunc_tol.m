function [U,S,V] = svdtrunc(A,tol)
% This code compute the fixed-accuracy lowrank factorization of the matrix A via the SVD
% of A.
%
% Copyright 2017 by Yingzhou Li and Haizhao Yang

[U,S,V] = svd(A,'econ');
if min(size(A)) > 0
    idx = find(diag(S)>tol*S(1,1)*max(size(A)));
    U = U(:,idx);
    S = S(idx,idx);
    V = V(:,idx);
end

end
