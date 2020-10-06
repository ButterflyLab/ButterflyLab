function [U,S,V] = svdtrunc(A,r,tol)

ms = min(size(A));
[U,S,V] = svd(A,'econ');
if(ms>r)
    idx = find(find(diag(S)>tol*S(1,1))<=r);
    U = U(:,idx);
    S = S(idx,idx);
    V = V(:,idx);
end

end
