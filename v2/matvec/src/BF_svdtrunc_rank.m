function [U,S,V] = BF_svdtrunc_rank(A,r,tol)
% compute the fixed rank approximation of A

ms = min(size(A));
[U,S,V] = svd(A,'econ');
if(ms>r)
    idx = find(find(diag(S)>tol*S(1,1))<=r);
    U = U(:,idx);
    S = S(idx,idx);
    V = V(:,idx);
end

end
