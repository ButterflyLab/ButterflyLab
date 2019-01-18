function [CU,sx,sx2] = CURL(fun,x,k,r,tol,opt)
% This code compute the CUR factorization using O(N) operations for the
% matrix fun(x,k) with an input rank estimation r and tolerence tol. "type"
% specify the output "sx".
%
% Output the CU factor and the selected row indices.
%
% Copyright
% 2018 modified by Qiyuan Pang and Haizhao Yang
% 2017 created by Yingzhou Li and Haizhao Yang

if nargin < 6, opt = 1; end

rr = 5*r;
% x = x(:);
% k = k(:);
xlen = size(x,1);
klen = size(k,1);

% idxk = randperm(klen);
if rr < klen
    grid = BF_Chey_grid(rr);
    idxk = round(grid*(klen-min(klen,rr)) + (0:min(klen,rr)-1)')+1;
    idxk= unique(idxk);
else
    idxk = 1:klen;
end
pk = k(idxk(1:min(klen,rr)),:);

Asub = fun(x,pk);

[~,Rmat,E] = qr(Asub',0);
if opt > 0
    if xlen*klen > 0
        Rmat = diag(Rmat);
        rr = find( abs(Rmat/Rmat(1)) > tol, 1, 'last');
        rr = min(r,rr);
    end
else
    rr = r;
end
sx2 = E(1:min(xlen,rr));
sx = x(sx2,:);
if 0
    [Uf,Sf,Vf] = BF_svdtrunc_tol(fun(sx,pk),tol);
    Sf = diag(Sf);
    CU = (fun(x,pk)*Vf)*diag(1./Sf)*Uf';
else
    CU = fun(x,pk)/fun(sx,pk);
end
end
