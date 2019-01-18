
%
% Interpolative Decomposition Butterfly Factorization, preprint, 2018.
% Copyright 2018 by Qiyuan Pang, Ken Ho, and Haizhao Yang

% Given submatrix factors S(1,1), S(2,1), S(1,2), S(2,2), form
%
%       [S(1,1)                     ]
%   S = [              S(2,1)       ]
%       [       S(1,2)              ]
%       [                     S(2,2)]
function S = BF_IDBF_merge_s(S2)
    assert(all(size(S2) == 2))  % check 2x2 size

    % work array
    e = cell(2);
    t = struct('I',e,'J',e,'V',e,'M',e,'N',e,'nz',e);

    % data for each splitting block
    for q = 1:2
    for p = 1:2
        [t(p,q).I,t(p,q).J,t(p,q).V] = find(S2{p,q});
        [t(p,q).M,t(p,q).N] = size(S2{p,q});
        t(p,q).nz = nnz(S2{p,q});
    end
    end
    nz = sum([t(:).nz]);  % total number of nonzeros

    % allocate sparse storage and fill
    I = zeros(nz,1);
    J = zeros(nz,1);
    V = zeros(nz,1);
    roff = [0 t(1,1).M+t(2,1).M; t(1,1).M t(1,1).M+t(2,1).M+t(1,2).M];
    coff = [0 t(1,1).N; t(1,1).N+t(1,2).N t(1,1).N+t(1,2).N+t(2,1).N];
    nz = 0;
    for q = 1:2
    for p = 1:2
        nv = t(p,q).nz;
        I(nz+(1:nv)) = roff(p,q) + t(p,q).I;
        J(nz+(1:nv)) = coff(p,q) + t(p,q).J;
        V(nz+(1:nv)) = t(p,q).V;
        nz = nz + nv;
    end
    end
    M = sum([t(:).M]);
    N = sum([t(:).N]);
    S = sparse(I,J,V,M,N);
end


