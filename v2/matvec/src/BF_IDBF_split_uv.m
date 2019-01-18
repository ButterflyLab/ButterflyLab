
%
% Interpolative Decomposition Butterfly Factorization, preprint, 2018.
% Copyright 2018 by Qiyuan Pang, Ken Ho, and Haizhao Yang

% Given a partition block
%
%   S(p,q) = [S11 S12] = [U11*S11*V11 U12*S12*V12]
%            [S21 S22]   [U21*S21*V21 U22*S22*V22]
%
% compute sparse U and V factors
%
%                                 [V11    ]
%   U = [U11     U12    ]  ,  V = [    V12]
%       [    U21     U22]         [V21    ]
%                                 [    V22]
%
% so that
%
%                [S11            ]
%   S(p,q) = U * [        S21    ] * V
%                [    S12        ]
%                [            S22]
function [U,V] = BF_IDBF_split_uv(data,p,q)
    [nb,np] = size(data);      % np is number of partitions after splitting
    nbp = nb/np;               % number of blocks per partition
    pp = 2*(p-1);              % pre-splitting partition index
    bp = pp*nbp+1:(pp+2)*nbp;  % block indices in row partition block

    % work array
    e = cell(2,1);
    t = struct('I',e,'J',e,'V',e,'M',e,'N',e);

    % row factors
    for i = 1:2          % loop over splitting index
        qq = 2*(q-1)+i;  % pre-splitting column partition index

        % count total nonzeros needed
        nz = 0;
        for b = bp                % loop over row blocks
            sk = data(b,qq).rsk;
            rd = data(b,qq).rrd;
            k  = length(sk);
            n  = k + length(rd);
            nz = nz + (n-k+1)*k;  % nonzeros to store [I T] for ID
        end

        % allocate sparse storage and fill U(:,i) = diag(U(1,i), ..., U(nb,i))
        I = zeros(nz,1);
        J = zeros(nz,1);
        V = zeros(nz,1);
        M  = 0;     % running row index
        N  = 0;     % running col index
        nz = 0;     % running number of nonzeros
        for b = bp  % loop over row blocks
            % extract ID interpolation factor
            sk = data(b,qq).rsk;
            rd = data(b,qq).rrd;
            T  = data(b,qq).rT;
            k = length(sk);
            n = k + length(rd);
            X = zeros(n,k);
            X(sk,:) = eye(k);
            X(rd,:) = T';
            % store in sparse matrix
            [I_,J_,V_] = find(X);
            nv = length(V_);
            I(nz+(1:nv)) = M + I_;
            J(nz+(1:nv)) = N + J_;
            V(nz+(1:nv)) = V_;
            M  = M  + n;
            N  = N  + k;
            nz = nz + nv;
        end
        % store for later access
        t(i).I = I;
        t(i).J = J;
        t(i).V = V;
        t(i).M = M;
        t(i).N = N;
    end

    % form 2x4 U factor U = [U(:,1) U(:,2)]
    N  = 0;
    nz = 0;
    for i = 1:2
        nv = length(t(i).V);
        I(nz+(1:nv)) = t(i).I;
        J(nz+(1:nv)) = t(i).J + N;
        V(nz+(1:nv)) = t(i).V;
        N  = N  + t(i).N;
        nz = nz + nv;
    end
    U = sparse(I,J,V,M,N);

    % column factors
    for i = 1:2          % loop over splitting index
        qq = 2*(q-1)+i;  % pre-splitting row partition index

        % count total nonzeros needed
        nz = 0;
        for b = bp                % loop over column blocks
            sk = data(b,qq).csk;
            rd = data(b,qq).crd;
            k  = length(sk);
            n  = k + length(rd);
            nz = nz + (n-k+1)*k;  % nonzeros to store [I T] for ID
        end

        % allocate sparse storage and fill V(:,i) = diag(V(1,i), ..., V(nb,1))
        I = zeros(nz,1);
        J = zeros(nz,1);
        V = zeros(nz,1);
        M  = 0;     % running row index
        N  = 0;     % running col index
        nz = 0;     % running number of nonzers
        for b = bp  % loop over column blocks
            % extract ID interpolation factor
            sk = data(b,qq).csk;
            rd = data(b,qq).crd;
            T  = data(b,qq).cT;
            k = length(sk);
            n = k + length(rd);
            X = zeros(k,n);
            X(:,sk) = eye(k);
            X(:,rd) = T;
            % store in sparse matrix
            [I_,J_,V_] = find(X);
            nv = length(V_);
            I(nz+(1:nv)) = M + I_;
            J(nz+(1:nv)) = N + J_;
            V(nz+(1:nv)) = V_;
            M  = M  + k;
            N  = N  + n;
            nz = nz + nv;
        end
        % store for later access
        t(i).I = I;
        t(i).J = J;
        t(i).V = V;
        t(i).M = M;
        t(i).N = N;
    end

    % form 4x2 V factor V = [V(:,1); V(:,2)]
    M  = 0;
    nz = 0;
    for i = 1:2
        nv = length(t(i).V);
        I(nz+(1:nv)) = t(i).I + M;
        J(nz+(1:nv)) = t(i).J;
        V(nz+(1:nv)) = t(i).V;
        M  = M  + t(i).M;
        nz = nz + nv;
    end
    V = sparse(I,J,V,M,N);
end

