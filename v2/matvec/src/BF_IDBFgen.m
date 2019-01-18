function [F] = BF_IDBFgen(fun_org,xx,kk,nn,rank_or_tol,flag)
% t, time for factorization
% ta, time for assembling
% F, the butterfly factorization of A
%
% Interpolative Decomposition Butterfly Factorization, preprint, 2018.
% Copyright 2018 by Qiyuan Pang, Ken Ho, and Haizhao Yang

if nargin < 6, flag = 0; end
%assert(mod(N2,1) == 0 && mod(n2,1) == 0)  % check integer input


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Dim = size(xx,2);
xbox = zeros(2,Dim);
kbox = zeros(2,Dim);
[Nx,~] = size(xx);
[Nk,~] = size(kk);
Nx = Nx^(1/Dim)*ones(1,Dim);
Nk = Nk^(1/Dim)*ones(1,Dim);
npx = 2.^round(log2(sqrt(Nx)));
npk = 2.^round(log2(sqrt(Nk)));
xbox(1,:) = min(xx);
xbox(2,:) = max(xx)+(max(xx)-min(xx))./Nx;
kbox(1,:) = min(kk);
kbox(2,:) = max(kk)+(max(kk)-min(kk))./Nk;
% make sure that the smallest dimension of the low-rank submatrices is
% larger than max(8,r)
levels = max(0,min(floor(log2([npx npk])-max(3,ceil(log2(nn))))));
Dim = size(xx,2);

xxboxidx = zeros(size(xx));
npxx = npx*2^levels;
xxset = cell(prod(npxx),1);
for i = 1:Dim
    edges = linspace(xbox(1,i),xbox(2,i),npxx(i)+1);
    [~,xxboxidx(:,i)] = histc(xx(:,i),edges);
end
[xxboxidx,xxidx] = sortrows(xxboxidx,Dim:-1:1);
[xxC,xxIA,~] = unique(xxboxidx,'rows','stable');
xxC = [xxC;zeros(1,size(xxC,2))];
xxIA = [xxIA;size(xx,1)+1];
for itx = 1:prod(npxx)
    x = BF_idx2vec(npxx,itx);
    if any(x ~= xxC(itx,:))
        xxC = [xxC(1:itx-1,:);x;xxC(itx:end,:)];
        xxIA = [xxIA(1:itx-1);xxIA(itx);xxIA(itx:end)];
    end
end
nx = size(xx,1);
iidx = zeros(nx,1);
ax = 0;
for itx = 1:prod(npxx)
    xxset{itx} = xxidx(xxIA(itx):xxIA(itx+1)-1)';
    for i = 1:length(xxset{itx})
        iidx(ax + i) = xxset{itx}(i);
    end
    ax = ax + length(xxset{itx});
end

kkboxidx = zeros(size(kk));
npkk = npk*2^levels;
kkset = cell(prod(npkk),1);
for i = 1:Dim
    edges = linspace(kbox(1,i),kbox(2,i),npkk(i)+1);
    [~,kkboxidx(:,i)] = histc(kk(:,i),edges);
end
[kkboxidx,kkidx] = sortrows(kkboxidx,Dim:-1:1);
[kkC,kkIA,~] = unique(kkboxidx,'rows','stable');
kkC = [kkC;zeros(1,size(kkC,2))];
kkIA = [kkIA;size(kk,1)+1];
for itk = 1:prod(npkk)
    k = BF_idx2vec(npkk,itk);
    if any(k ~= kkC(itk,:))
        kkC = [kkC(1:itk-1,:);k;kkC(itk:end,:)];
        kkIA = [kkIA(1:itk-1);kkIA(itk);kkIA(itk:end)];
    end
end
nk = size(kk,1);
iidk = zeros(nk,1);
ak = 0;
for itk = 1:prod(npkk)
    kkset{itk} = kkidx(kkIA(itk):kkIA(itk+1)-1)';
    for i = 1:length(kkset{itk})
        iidk(ak + i) = kkset{itk}(i);
    end
    ak = ak + length(kkset{itk});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A is a function handle for evaluating A(i,j) via A(i,j)
%N  = 2^N2;              % total size
%n  = 2^n2;              % leaf block size
%L  = floor((N2-n2)/2);  % number of compression levels
%nb = N/n;               % number of leaf blocks (= # to compress at each level)
fun = @(i,j)fun_org(xx(i,:),kk(j,:));
N = size(xx,1);
nb = prod(npxx);
L = floor(log2(nb)/2);

% print input summary
if flag, fprintf('total size = %d, leaf size = %d, levels = %d',N,n,L); end
if rank_or_tol < 1, s = 'tol'; else s = 'rank'; end
if flag, fprintf(', %s = %g\n',s,rank_or_tol); end

% Matrix to compress: just do the standard Fourier matrix for simplicity. By
% symmetry, we can use the same uniform binary tree on both the spatial and
% frequency domains. The matrix is specified through a function to generate
% arbitrary submatrices in order to allow efficient subsampling. In a more
% complete code, we would probably need a "proxy" kernel function as well.
% A = @(i,j)(exp(-2i*pi*(i(:)-1)*(j(:)'-1)/N));

F.data = cell(L,1);  % set up storage for compression data

% We organize the complementary butterfly structure in terms of "blocks" and
% "partitions" at each compression level, where the former indexes one domain at
% the lower tree level (closer to the leaves) and the latter the other domain at
% the complementary higher level (closer to the root). Each block will be
% compressed with respect to each partition hierarchically; this applies
% independently to both the rows and columns. The compression will be done using
% the ID, so we need to keep track of which row and column indices remain in
% each block for each partition; these will get thinned out as we move up the
% tree. The variables nb and np will track the numbers of blocks and partitions,
% respectively.
np = 1;
ridx = cell(nb,np);
cidx = cell(nb,np);
for b = 1:nb
    ridx{b,1} = xxset{b};  % initial indices in each block
    cidx{b,1} = kkset{b};        % same tree structure as rows
end


nb_ = nb;  % save initial number of blocks (= # total over all partitions)

% Do the compression. Our basic approach here essentially just generates the
% whole matrix and so has O(N^2) complexity. A more general code could bring
% this down to O(N*log(N)) by using proxy points (e.g., Chebyshev), but we did
% not do this here for simplicity. However, in this particular example, the
% matrix is translation-invariant within each partition (actually across all
% partitions, up to twiddle factors), so we can "cheat" by compressing only one
% block per partition. This recovers O(N*log(N)) complexity.
cheat = 0;  % "cheat" by using translation invariance?
if flag, fprintf('---- start compression ----\n'); end
%start = tic;
for l = 1:L  % loop from leaf level up
%    tic;
    
    % aggregate for next level
    ridx_ = ridx;          % work arrays
    cidx_ = cidx;
    nb = nb/2;             % butterfly split: halve number of blocks and double
    ridx = cell(nb,2*np);  %   number of partitions; equivalent to doubling and
    cidx = cell(nb,2*np);  %   halving respective domains
    for p = 1:np           % loop over old number of partitions
        for b = 1:nb       % loop over new number of blocks
            % update remaining indices by pulling from appropriate data
            ridx{b,2*p-1} = [ridx_{2*b-1,p} ridx_{2*b,p}];
            cidx{b,2*p-1} = [cidx_{2*b-1,p} cidx_{2*b,p}];
            ridx{b,2*p} = ridx{b,2*p-1};
            cidx{b,2*p} = cidx{b,2*p-1};
        end
    end
    np = np*2;  % update number of partitions (postponed for easier indexing)
    
    % set up storage for each level
    e = cell(nb,np);
    F.data{l} = struct('ri',e,'rsk',e,'rrd',e,'rT',e,...   % row data
        'ci',e,'csk',e,'crd',e,'cT',e);  % col data
    
    % number of blocks at current level making up each partition; used to get
    % "complementary" indices by reference to remaining indices at current level
    nbp = nb/np;
    
    % row compression
    for p = 1:np                                 % loop over partitions
        ci = [cidx{(p-1)*nbp+1:p*nbp,p}];        % complementary cols
        for b = 1:nb                             % loop over blocks
            ri = ridx{b,p};                      % rows to compress
            if ~cheat || b == 1
                K = fun(ri,ci);                    % pull out submatrix ...
                [sk,rd,T] = BF_ID(K',rank_or_tol);  % ... and compress rows with ID
            end
            F.data{l}(b,p).ri  = ri;             % store row block indices
            F.data{l}(b,p).rsk = sk;             % store skeleton indices
            F.data{l}(b,p).rrd = rd;             % store redundant indices
            F.data{l}(b,p).rT  = T;              % store interpolation matrix
            ridx{b,p} = ri(sk);                  % restrict to skeletons
        end
    end
    
    % column compression
    for p = 1:np                                 % loop over partitions
        ri = [ridx{(p-1)*nbp+1:p*nbp,p}];        % complementary rows
        for b = 1:nb                             % loop over blocks
            ci = cidx{b,p};                      % cols to compress
            if ~cheat || b == 1
                K = fun(ri,ci);                    % pull out submatrix ...
                [sk,rd,T] = BF_ID(K,rank_or_tol);   % ... and compress cols with ID
            end
            F.data{l}(b,p).ci  = ci;             % store col block indices
            F.data{l}(b,p).csk = sk;             % store skeleton indices
            F.data{l}(b,p).crd = rd;             % store redundant indices
            F.data{l}(b,p).cT  = T;              % store interpolation matrix
            cidx{b,p} = ci(sk);                  % restrict to skeletons
        end
    end
    
    % print compression statistics: average block size -> average # skeletons
%    t = toc;
%    cn = (length([F.data{l}(:).ri ]) + length([F.data{l}(:).ci ]))/(2*nb_);
%    ck = (length([F.data{l}(:).rsk]) + length([F.data{l}(:).csk]))/(2*nb_);
%    if flag, fprintf('level %2d: %6.2f x %8.2f -> %6.2f in %8.4e s\n',l,cn,cn*nbp,ck,t); end
end
%if flag, fprintf('---- done compression: %8.4e s ----\n',toc(start)); end

%start = tic;  % timer for sparse matrix assembly

% Reform into sparse matrix factors: first, the left and right matrices U and V.
% Recall that the sparse matrix factorization is based on the recursive
% splitting of a given "partition blocK"
%
%   S = [S11 S12] = [U11*S11*V11 U12*S12*V12]
%       [S21 S22]   [U21*S21*V21 U22*S22*V22]
%
% as
%
%                     [S11            ] [V11    ]
%   [U11     U12    ] [        S21    ] [    V12]
%   [    U21     U22] [    S12        ] [V21    ]
%                     [            S22] [    V22]
%
% This part of the code does two things at each level: (1) compute the matrices
% U(p,q) and V(p,q) as above to split each partition block S(p,q) into 4x4
% blocks; and (2) diagonally merge the U(p,q) and V(p,q) to form global matrices
% U and V for the whole level.
%
% Note that the notation for V below is "backwards": to keep the same logic for
% U and V, we always let the first index range over the relevant dimension
% (i.e., rows or columns, as appropriate) and the second index over the
% complementary other dimension. In other words, U12 in the example above might
% be U{1,2}, while V12 is V{2,1}.
F.U = cell(L,1);
F.V = cell(L,1);

e = cell(2,1);
t = struct('I',e,'J',e,'V',e,'totalM',e,'totalN',e);

for i = 1:size(F.data{1},2) 
    nz = 0;
    for j = 1:size(F.data{1},1)
        rd = F.data{1}(j,i).rrd;
        sk = F.data{1}(j,i).rsk;
        k = length(sk);
        n = k + length(rd);
        nz = nz + (n-k+1)*k;
    end
    
    I = zeros(nz,1);
    J = zeros(nz,1);
    V = zeros(nz,1);
    totalM = 0;
    totalN = 0;
    nz = 0;
    for j = 1:size(F.data{1},1) 
        rd = F.data{1}(j,i).rrd;
        sk = F.data{1}(j,i).rsk;
        T = F.data{1}(j,i).rT;
        k = length(sk);
        n = k + length(rd);
        X = zeros(n,k);
        X(sk,:) = eye(k);
        X(rd,:) = T';
        [I_,J_,V_] = find(X);
        nv = length(V_);
        I(nz+(1:nv)) = iidx(totalM + I_);
        J(nz+(1:nv)) = totalN + J_;
        V(nz+(1:nv)) = V_;
        totalM  = totalM  + n;
        totalN  = totalN  + k;
        nz = nz + nv;
    end
    t(i).I = I;
    t(i).J = J;
    t(i).V = V;
    t(i).totalM = totalM;
    t(i).totalN = totalN;
end
NN  = 0;
nz = 0;
for i = 1:2
    nv = length(t(i).V);
    I(nz+(1:nv)) = t(i).I;
    J(nz+(1:nv)) = t(i).J + NN;
    V(nz+(1:nv)) = t(i).V;
    NN  = NN  + t(i).totalN;
    nz = nz + nv;
end
F.U{1} = sparse(I,J,V,totalM,NN);

for i = 1:size(F.data{1},2) 
    nz = 0;
    for j = 1:size(F.data{1},1)
        rd = F.data{1}(j,i).crd;
        sk = F.data{1}(j,i).csk;
        k = length(sk);
        n = k + length(rd);
        nz = nz + (n-k+1)*k;
    end
    
    I = zeros(nz,1);
    J = zeros(nz,1);
    V = zeros(nz,1);
    totalM = 0;
    totalN = 0;
    nz = 0;
    for j = 1:size(F.data{1},1) 
        rd = F.data{1}(j,i).crd;
        sk = F.data{1}(j,i).csk;
        T = F.data{1}(j,i).cT;
        k = length(sk);
        n = k + length(rd);
        X = zeros(k,n);
        X(:,sk) = eye(k);
        X(:,rd) = T;
        [I_,J_,V_] = find(X);
        nv = length(V_);
        I(nz+(1:nv)) = totalM + I_;
        J(nz+(1:nv)) = iidk(totalN + J_);
        V(nz+(1:nv)) = V_;
        totalM  = totalM  + k;
        totalN  = totalN  + n;
        nz = nz + nv;
    end
    t(i).I = I;
    t(i).J = J;
    t(i).V = V;
    t(i).totalM = totalM;
    t(i).totalN = totalN;
end
MM  = 0;
nz = 0;
for i = 1:2
    nv = length(t(i).V);
    I(nz+(1:nv)) = t(i).I+MM;
    J(nz+(1:nv)) = t(i).J;
    V(nz+(1:nv)) = t(i).V;
    MM  = MM  + t(i).totalM;
    nz = nz + nv;
end
F.V{1} = sparse(I,J,V,MM,totalN);


for l = 2:L  % loop over levels
    
    % compute U and V factors for each partition block
    np = size(F.data{l},2)/2;
    U = cell(np,np);%%%%%cell(np)?
    V = cell(np,np);%%%%%cell(np)?
    for q = 1:np  % loop over all partition blocks ...
        for p = 1:np  % ... and compute U and V splitting factors
            [U{p,q},V{p,q}] = BF_IDBF_split_uv(F.data{l},p,q);
        end
    end
    
    % recursively merge U and V factors over all partition blocks
    while np > 1
        np = np/2;
        for q = 1:np  % loop over all partition blocks ...
            for p = 1:np  % ... and do 2x2 merge until no more
                U{p,q} = BF_IDBF_merge_uv(U(2*p-1:2*p,2*q-1:2*q));
                V{p,q} = BF_IDBF_merge_uv(V(2*p-1:2*p,2*q-1:2*q));
            end
        end
        U = U(1:np,1:np);
        V = V(1:np,1:np);
    end
    
    % store merged U and V factors
    F.U{l} = U{1};
    F.V{l} = V{1};
end

% Now form middle sparse S factor. This takes the form of the 4x4 block
% splitting above, where each block is then recursively split in the same way.
% First, just form the submatrix of S corresponding to each partition block.
[nb,np] = size(F.data{L});
nbp = nb/np;
S   = cell(np);
tmp = cell(nbp,1);  % work array
for q = 1:np                     % loop over secondary partition blocks
    bq = (q-1)*nbp+1:q*nbp;      % blocks in secondary partition
    for p = 1:np                 % loop over primary partition blocks
        bp = (p-1)*nbp+1:p*nbp;  % blocks in primary partition
        for b = 1:nbp
            ri = F.data{L}(bp(b),q).ri;
            sk = F.data{L}(bp(b),q).rsk;
            tmp{b} = ri(sk);
        end
        ri = [tmp{:}];           % all remaining rows in (p,q) block
        for b = 1:nbp
            ci = F.data{L}(bq(b),p).ci;
            sk = F.data{L}(bq(b),p).csk;
            tmp{b} = ci(sk);
        end
        ci = [tmp{:}];           % all remaining cols in (p,q) block
        S{p,q} = fun(ri,ci);       % generate skeleton submatrix
    end
end

% Then recursively merge S factors until we have one global S matrix.
for l = 1:L
    np = np/2;
    for q = 1:np
        for p = 1:np
            S{p,q} = BF_IDBF_merge_s(S(2*p-1:2*p,2*q-1:2*q));
        end
    end
    S = S(1:np,1:np);
end
F.S = S{1};  % store factorization

% print sparse matrix statistics
%ta = toc(start);
nz = nnz(F.S);
for l = 1:L
    nz = nz + nnz(F.U{l}) + nnz(F.V{l});
end
if flag, 
fprintf('sparse matrix assembly in %8.4e s, nnz = %d\n',ta,nz); end
end