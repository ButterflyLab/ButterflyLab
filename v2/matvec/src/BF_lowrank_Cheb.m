function [U,V] = BF_lowrank_Cheb(x,k,fun,type,tol,mR,isSVD,tR)
% Input
% x -     samples in x
% k -     samples in k
% fun -   a function handle in (x,k)
% type -  1: interpolation in x
%         2: interpolation in k
% tol -   low-rank approximation accuracy
% tR -    test rank of the low-rank approximation
% mR -    maximum rank of the low-rank approximation
% isSVD - whether use SVD to compress the low-rank factorization by
%         Chebyshev interpolation
%
% Output
% U*V' \approx fun(x,k)
%
% Copyright 2018 by Haizhao Yang, National University of Singapore

if nargin < 7, isSVD = 1; end
if nargin < 6, mR = log2(numel(x)); end
if nargin < 5, tol = 1e-13; end;
if nargin < 4, type = 2; end;
if nargin < 8, tR = mR; end

grid = BF_Chey_grid(tR);
switch type
    case 1 % interpolation in x
        Mx = max(x); mx = min(x);
        grid = grid*(Mx-mx) + mx;
        U = BF_Lagrange(grid(:),x(:));
        U = U';
        V = fun(grid(:),k(:));
        if isSVD
            [V,S,U2] = svd(V','econ');
            if tol > 0
                pos = find(diag(S(1:mR,1:mR))>S(1,1)*tol);
            else
                pos = 1:mR;
            end
            U = U*U2(:,pos)*S(pos,pos);
            V = V(:,pos);
        else
            V = V';
        end
    case 2 % interpolation in k
        Mk = max(k); mk = min(k);
        grid = grid*(Mk-mk) + mk;
        V = BF_Lagrange(grid(:),k(:));
        U = fun(x(:),grid(:));
        if isSVD
            [U,S,V2] = svd(U,'econ');
            if tol > 0
                pos = find(diag(S(1:mR,1:mR))>S(1,1)*tol);
            else
                pos = 1:mR;
            end
            U = U(:,pos);
            V = V'*V2(:,pos)*S(pos,pos);
        else
            V = V';
        end
end
end
