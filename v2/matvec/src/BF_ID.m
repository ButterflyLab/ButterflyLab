% Compute ID approximation A(:,rd) ~ A(:,sk)*T. The precision is specified by
% rank_or_tol; if rank_or_tol < 1, then it is taken as the relative precision
% and otherwise as the rank.
%
% This code is slow because we don't use subsampling techniques
%
% Copyright 2018 by Qiyuan Pang, Ken Ho, and Haizhao Yang

function [sk,rd,T] = BF_ID(A,rank_or_tol)
    [m,n] = size(A);
    if isempty(A)
        sk = [];
        rd = 1:n;
        T = zeros(0,n);
        return
    end
    [~,R,E] = qr(A,0);
    if rank_or_tol < 1
        k = sum(abs(diag(R)) > abs(R(1))*rank_or_tol);
    else
        k = min(rank_or_tol,n);
    end
    sk = E(1:k);
    rd = E(k+1:end);
    T = R(1:k,1:k)\R(1:k,k+1:end);
end