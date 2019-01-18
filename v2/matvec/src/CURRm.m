function [UR,sk] = CURR(fun,x,k,r,tol,type)
% This code compute the CUR factorization using O(N) operations for the
% matrix fun(x,k) with an input rank estimation r and tolerence tol. "type"
% specify the output "sx".
%
% Output the UR factor and the selected column indices.
%
% Copyright
% 2018 modified by Qiyuan Pang and Haizhao Yang
% 2017 created by Yingzhou Li and Haizhao Yang

if nargin < 6, type = 1; end

rr = 5*r;
% x = x(:);
% k = k(:);
xlen = size(x,1);
klen = size(k,1);

idxx = randperm(xlen);
px = x(idxx(1:min(xlen,rr)),:);

Asub = fun(px,k);

[~,Rmat,E] = qr(Asub,0);
if xlen*klen > 0
    Rmat = diag(Rmat);
    rr = find( abs(Rmat/Rmat(1)) > tol, 1, 'last');
    rr = min(r,rr);
end
sk2 = E(1:min(klen,rr));
sk = k(sk2,:);
if 0
    [Uf,Sf,Vf] = BF_svdtrunc_tol(fun(px,sk),tol);
    Sf = diag(Sf);
    UR = Vf*diag(1./Sf)*(Uf'*fun(px,k));
else
    UR = fun(px,sk)\fun(px,k);
end
if type == 2
    sk = sk2;
end
end
