function [U,S,V] = BF_rsvd(L,R,K)
% randomized SVD for the low-rank matrix L*R

M = size(L,1);
N = size(R,2);
P = min(2*K,N);
X = randn(N,P);
Y = L*(R*X);
W1 = orth(Y);
B = (W1'*L)*R;
[W2,S,V] = svd(B,'econ');
U = W1*W2;
K=min(K,size(U,2));
U = U(:,1:K);
S = S(1:K,1:K);
V=V(:,1:K);