function X = HSSBF_Schulz_mat(A,X,maxIter)
% Input:
% A: the matrix to be inverted
% X: an approximate inverse of A
% maxIter: the maximumn iteration number
%
% Output:
% X: a better approximate inverse of A using the Schultz iteration

if nargin < 2, maxIter = 5; end
for cnt = 1:maxIter
    X = 2*X-X*A*X;
end
end