function [U,sx,sk] = CUR(fun,x,k,r,tol)
% This code compute the CUR factorization using O(1) operations for the
% matrix fun(x,k) with an input rank estimation r.
%
% Only the middle factor U and the selected row/column indices are output.
%
% Copyright 2017 by Yingzhou Li and Haizhao Yang

rr = 5*r;
% x = x(:);
% k = k(:);
xlen = size(x,1);
klen = size(k,1);

idxx = randperm(xlen);
idxk = randperm(klen);
px = x(idxx(1:min(xlen,rr)),:);
pk = k(idxk(1:min(klen,rr)),:);

Asub = fun(px,pk);

[~,Rmat,E] = qr(Asub,0);
if xlen*klen > 0
    Rmat = diag(Rmat);
    rr = find( abs(Rmat/Rmat(1)) > tol, 1, 'last');
    rr = min(r,rr);
end
sk = pk(E(1:min(klen,rr)),:);

[~,Rmat,E] = qr(Asub',0);
if xlen*klen > 0
    Rmat = diag(Rmat);
    rr = find( abs(Rmat/Rmat(1)) > tol, 1, 'last');
    rr = min(r,rr);
end
sx = px(E(1:min(xlen,rr)),:);

[Uf,Sf,Vf] = BF_svdtrunc_tol(fun(sx,sk),tol);
Sf = diag(Sf);
U = Vf*diag(1./Sf)*Uf';

end
