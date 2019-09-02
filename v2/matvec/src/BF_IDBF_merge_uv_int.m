% Multidimensional Phase Recovery and Interpolative Decomposition 
% Butterfly Factorization, preprint, 2019.
% Copyright 2018 by Ze Chen, Juan Zhang, Ken Ho, and Haizhao Yang

% Given splitting factors U(1,1), ..., U(2^Dim,2^Dim) for 
% S(1,1), ..., S(2^Dim,2^Dim), form
%
%       [U_1                  ]           [U(1,k)                      ]
%   U = [    U_2              ] with Uk = [       U(2,k)               ]
%       [        ...          ]           [              ...           ]
%       [            U_{2^Dim}]           [                  U(2^Dim,k)]
%
% See recursive splitting formula.
% The same code works for V(1,1), ..., V(2^Dim,2^Dim).
% provided that the first index is over the columns (primary dimension) 
% and the second index over the rows (complementary secondary dimension).

function U = BF_IDBF_merge_uv_int(U2,Dim)
    assert(size(U2,1) == (2^Dim) & size(U2,2) == (2^Dim))
    % check 2^Dim * 2^Dim size

    % work array
    e = cell(2^Dim);
    t = struct('I',e,'J',e,'V',e,'M',e,'N',e,'nz',e);

    % data for each splitting block
    for q = 1:(2^Dim)
        for p = 1:(2^Dim)
            [t(p,q).I,t(p,q).J,t(p,q).V] = find(U2{p,q});
            [t(p,q).M,t(p,q).N] = size(U2{p,q});
            t(p,q).nz = length(t(p,q).V);
        end
    end
    nz = sum([t(:).nz]);  % total number of nonzeros

    % allocate sparse storage and fill
    I = zeros(nz,1);
    J = zeros(nz,1);
    V = zeros(nz,1);
    M  = 0;
    N  = 0;
    nz = 0;
    for q = 1:(2^Dim)
    for p = 1:(2^Dim)
        nv = t(p,q).nz;
        I(nz+(1:nv)) = M + t(p,q).I;
        J(nz+(1:nv)) = N + t(p,q).J;
        V(nz+(1:nv)) =     t(p,q).V;
        M  = M  + t(p,q).M;
        N  = N  + t(p,q).N;
        nz = nz + nv;
    end
    end
    U = sparse(I,J,V,M,N);
end
