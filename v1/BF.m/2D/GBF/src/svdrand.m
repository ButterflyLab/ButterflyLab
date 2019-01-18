function [U,S,V] = svdrand(A,r,tol)

ms = min(size(A));
[m,n] = size(A);

if(ms>r)
    tR = r+3;
    Rcol = randn(n,tR);
    Rrow = randn(m,tR);
    Rcol = A*Rcol;
    Rrow = A'*Rrow;
    for iter = 1:0
        Rcol = A*(A'*Rcol);
        Rrow = A'*(A*Rrow);
    end
    [Qcol,~] = qr(Rcol,0);
    [Qrow,~] = qr(Rrow,0);
    [U,S,V] = svd(Qcol'*A*Qrow,'econ');
    idx = find(find(diag(S)>tol*S(1,1))<=r);
    U = Qcol*U(:,idx);
    S = S(idx,idx);
    V = Qrow*V(:,idx);
else
    [U,S,V] = svd(A,'econ');
    if(ms>0)
        idx = find(find(diag(S)>tol*S(1,1))<=r);
        U = U(:,idx);
        S = S(idx,idx);
        V = V(:,idx);
    end
end

end