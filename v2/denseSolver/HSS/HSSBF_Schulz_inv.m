function invFunNew = HSSBF_Schulz_inv(Factor,invFunOld,maxIter)
% Input:
% Factor: a structure storing the data-sparse representation of the inverse
%       of A(x,k), and some auxiliary variables.
% Factor.isLeaf: 1, if all matrices are dense
%                0, if all matrices are data-sparse
% Factor.sz: the size of A
% Factor.szsub: the size of A11
% Factor.B
% Factor.Btr
% Factor.Bctr
% Factor.A11
% Factor.A22
% Factor.A11inv
% Factor.A22inv
%
% invFunOld: a function handle representing the inverse of a HSSBF matrix
% using Schultz iteration
% maxIter - the maximum iteration number of the Schultz iteration
%
% Output:
% invFunNew: a function handle representing the inverse of a HSSBF matrix
% using Schultz iteration

if nargin < 3, maxIter = 5; end
if maxIter == 0
    invFunNew = @(f) HSSBF_sol(Factor,f);
else
    invFun = HSSBF_Schulz_inv(Factor,invFunOld,maxIter-1);
    invFunNew = @(f) applyFun(invFun,Factor,f);
end
end

function y = applyFun(invFunOld,Factor,f)
g = 2*invFunOld(f);
y = g-invFunOld(HSSBF_sol(Factor,g));
end