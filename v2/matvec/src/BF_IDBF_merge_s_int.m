% Multidimensional Phase Recovery and Interpolative Decomposition 
% Butterfly Factorization, preprint, 2019.
% Copyright 2018 by Ze Chen, Juan Zhang, Ken Ho, and Haizhao Yang

% Given submatrix factors S(1,1), ..., S(2^Dim,2^Dim), form
%
%       [  S_{1,1}     S_{1,2}   ...   S_{1,2^Dim}  ]
%   S = [  S_{2,1}     S_{2,2}   ...   S_{2,2^Dim}  ]
%       [    ...         ...     ...       ...      ]
%       [S_{2^Dim,1} S_{2^Dim,2} ... S_{2^Dim,2^Dim}]
% 
% with S_{i,j} as a 2^Dim * 2^Dim block matrix with the (j,i)-th block
% as S(j,i).

function S = BF_IDBF_merge_s_int(S2,Dim)
    assert(all(size(S2) == 2^Dim))  % check 2^Dim * 2^Dim size

    % work array
    e = cell(2^Dim);
    t = struct('I',e,'J',e,'V',e,'M',e,'N',e,'nz',e);

    % data for each splitting block
    M = zeros(2^Dim); N = zeros(2^Dim);
    for q = 1:(2^Dim)
        for p = 1:(2^Dim)
            [t(p,q).I,t(p,q).J,t(p,q).V] = find(S2{p,q});
            [t(p,q).M,t(p,q).N] = size(S2{p,q});
            M(p,q) = t(p,q).M; N(p,q) = t(p,q).N;
            t(p,q).nz = nnz(S2{p,q});
        end
    end
    nz = sum([t(:).nz]);  % total number of nonzeros

    % allocate sparse storage and fill
    I = zeros(nz,1);
    J = zeros(nz,1);
    V = zeros(nz,1);
    roff = zeros(2^Dim);
    for i = 2:(2^Dim^2)
        roff(i) = M(i-1)+roff(i-1);
    end
    coff = zeros(2^Dim); N = N';
    for i = 2:(2^Dim^2)
        coff(i) = N(i-1)+coff(i-1);
    end
    coff = coff';
    nz = 0;
    for q = 1:(2^Dim)
    for p = 1:(2^Dim)
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
