
%
% Interpolative Decomposition Butterfly Factorization, preprint, 2018.
% Copyright 2018 by Qiyuan Pang, Ken Ho, and Haizhao Yang

% Given splitting factors U(1,1), ..., U(2,2) for S(1,1), ..., S(2,2), form
%
%       [U(1,1)                     ]
%   U = [       U(2,1)              ]
%       [              U(1,2)       ]
%       [                     U(2,2)]
%
% See recursive splitting formula. The same code works for V(1,1), ..., V(2,2)
% provided that the first index is over the columns (primary dimension) and the
% second index over the rows (complementary secondary dimension).
function U = BF_IDBF_merge_uv(U2)
    assert(size(U2,1) == 2 & size(U2,2) == 2)  % check 2x2 size

    % work array
    e = cell(2);
    t = struct('I',e,'J',e,'V',e,'M',e,'N',e,'nz',e);

    % data for each splitting block
    for q = 1:2
    for p = 1:2
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
    for q = 1:2
    for p = 1:2
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


