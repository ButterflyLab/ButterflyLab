% Multidimensional Phase Recovery and Interpolative Decomposition 
% Butterfly Factorization, preprint, 2019.
% Copyright 2018 by Ze Chen, Juan Zhang, Ken Ho, and Haizhao Yang

% O(N log N) operation and memory complexity.

function P = BF_RecoveryPath(X)
% X - locations of the grid points
% P - recovery path matrix

    DT = delaunayTriangulation(X);
    
    E = edges(DT);
    d = sqrt(sum((X(E(:,1), :) - X(E(:,2), :)).^2, 2));
    G = graph(E(:,1), E(:,2), d);
    
    [Tree, ~] = minspantree(G);
    P = bfsearch(Tree, 1, 'edgetonew');
    
end