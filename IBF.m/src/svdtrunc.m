function [U,S,V] = svdtrunc(A,tol)

[U,S,V] = svd(A,'econ');
if min(size(A)) > 0
    idx = find(diag(S)>tol*S(1,1)*max(size(A)));
    U = U(:,idx);
    S = S(idx,idx);
    V = V(:,idx);
end

end
