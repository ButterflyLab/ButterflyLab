function F = BF_MIDBF_nonuniform(fun,xx,kk,nn,rk,tol,rand_or_cheb,tt,opt)
% fun - function handle such that fun(xx,kk) is the complementary
% low-rank matrix A to be compressed
% xx - column vector and must be (1:numel(xx))'
% kk - column vector and must be (1:numel(kk))'
% nn - smalles number of points in leaves
% rk - rank parameter
% tol - set up accuracy
% tt - over sampling parameter in ID
% opt - whether use adaptive rank or fix rank
%       if opt = 0, fix rank rk in the BF no matter what tol is
%       if opt = 1, use a rank less than or equal to rk trying to obtain
%        an accuracy as good as tol; may not be able to achieve tol if rk
%        is too small
% rand_or_cheb - whether use random sampling or Mock-Chebshev points for ID
% F, the butterfly factorization of A
%
%
% Multidimensional Phase Recovery and Interpolative Decomposition
% Butterfly Factorization, preprint, 2019.
% Copyright 2018 by Ze Chen, Juan Zhang, Ken Ho, and Haizhao Yang


if nargin < 7, rand_or_cheb = 'rand'; end
if nargin < 8, tt = 5; end
if nargin < 9, opt = 1; end % 1: adaptive rank; 0: exact rank=rk

if rk > 0
    isRk = 1; % use rk in the low-rank approximation
    grid = BF_Chey_grid(rk);
else
    isRk = 0;
end
if tol < 1
    isTol = 1; % use tol in the low-rank approximation
else
    isTol = 0;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
xorder = BF_Zorder(prod(npxx), Dim);
for itx = 1:prod(npxx)
    xxset{itx} = xxidx(xxIA(xorder(itx)):xxIA(xorder(itx)+1)-1)';
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
korder = BF_Zorder(prod(npkk), Dim);
for itk = 1:prod(npkk)
    kkset{itk} = kkidx(kkIA(korder(itk)):kkIA(korder(itk)+1)-1)';
    for i = 1:length(kkset{itk})
        iidk(ak + i) = kkset{itk}(i);
    end
    ak = ak + length(kkset{itk});
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = size(xx,1);
nb = prod(npxx);
L = floor(log2(nb)/(2*Dim));

% Matrix to compress

F.data = cell(L,1);  % set up storage for compression data

np = 1;
ridx = cell(nb,np);
cidx = cell(nb,np);
for b = 1:nb
    ridx{b,1} = xxset{b};  % initial indices in each block
    cidx{b,1} = kkset{b};  % same tree structure as rows
end

cheat = 0;  % "cheat" by using translation invariance?
for l = 1:L  % loop from leaf level up
    % aggregate for next level
    ridx_ = ridx;          % work arrays
    cidx_ = cidx;
    nb = nb/(2^Dim);  % butterfly split: halve number of blocks and double
    ridx = cell(nb,(2^Dim)*np);  % number of partitions; equivalent to
    cidx = cell(nb,(2^Dim)*np);  % doubling and halving respective domains
    for p = 1:np  % loop over old number of partitions
        for b = 1:nb  % loop over new number of blocks
            % update remaining indices by pulling from appropriate data
            rridx = []; ccidx = [];
            for blk = 1:(2^Dim)
                rridx = [rridx ridx_{(2^Dim)*(b-1)+blk,p}];
                ccidx = [ccidx cidx_{(2^Dim)*(b-1)+blk,p}];
            end
            for blk = 1:(2^Dim)
                ridx{b,(2^Dim)*(p-1)+blk} = rridx;
                cidx{b,(2^Dim)*(p-1)+blk} = ccidx;
            end
        end
    end
    np = np*(2^Dim);  % update number of partitions ...
    % ... (postponed for easier indexing)
    
    % set up storage for each level
    e = cell(nb,np);
    F.data{l} = struct('ri',e,'rsk',e,'rrd',e,'rT',e,...   % row data
        'ci',e,'csk',e,'crd',e,'cT',e);  % col data
    
    % number of blocks at current level making up each partition;
    % used to get "complementary" indices by reference to remaining
    % indices at current level
    nbp = nb/np;
    
    % row compression
    for p = 1:np  % loop over partitions
        for b1 = 1:np
            ci = [cidx{(p-1)*nbp+1:p*nbp,b1}];  % complementary cols
            for b2 = 1:nbp                    % loop over blocks
                b = (b1-1)*nbp + b2;
                ri = ridx{b,p};  % rows to compress
                if ~cheat || b == 1
                    if isRk && isTol
                        switch rand_or_cheb
                            case 'cheb'
                                [T,sk,~,rd] = BF_ID_Cheby(fun,ri(:),...
                                    ci(:),grid,rk,tol,'r',opt);
                            case 'rand'
                                [T,sk,~,rd] = BF_ID_rand(fun,ri(:),...
                                    ci(:),rk,tol,'r',tt,opt);
                        end
                    elseif isRk
                        switch rand_or_cheb
                            case 'cheb'
                                [T,sk,~,rd] = BF_ID_Cheby(fun,ri(:),...
                                    ci(:),grid,rk,1e-15,'r',opt);
                            case 'rand'
                                [T,sk,~,rd] = BF_ID_rand(fun,ri(:),...
                                    ci(:),rk,1e-15,'r',tt,opt);
                        end
                    else
                        switch rand_or_cheb
                            case 'cheb'
                                grid = BF_Chey_grid(size(ri(:),1));
                                [T,sk,~,rd] = BF_ID_Cheby(fun,ri(:),...
                                    ci(:),grid,size(ri(:),1),tol,'r',opt);
                            case 'rand'
                                [T,sk,~,rd] = BF_ID_rand(fun,ri(:),...
                                    ci(:),size(ri(:),1),tol,'r',tt,opt);
                        end
                    end
                end
                F.data{l}(b,p).ri  = ri;  % store row block indices
                F.data{l}(b,p).rsk = sk;  % store skeleton indices
                F.data{l}(b,p).rrd = rd;  % store redundant indices
                T(T==0) = eps;
                F.data{l}(b,p).rT  = T;   % store interpolation matrix
                ridx{b,p} = ri(sk);  % restrict to skeletons
            end
        end
    end
    
    % column compression
    for p = 1:np  % loop over partitions
        for b1 = 1:np
            ri = [ridx{(p-1)*nbp+1:p*nbp,b1}];    % complementary rows
            for b2 = 1:nbp                    % loop over blocks       
                b = (b1-1)*nbp + b2;
                ci = cidx{b,p};  % cols to compress
                if ~cheat || b == 1
                    if isRk && isTol
                        switch rand_or_cheb
                            case 'cheb'
                                [T,sk,~,rd] = BF_ID_Cheby(fun,ri(:),...
                                    ci(:),grid,rk,tol,'c',opt);
                            case 'rand'
                                [T,sk,~,rd] = BF_ID_rand(fun,ri(:),...
                                    ci(:),rk,tol,'c',tt,opt);
                        end
                    elseif isRk
                        switch rand_or_cheb
                            case 'cheb'
                                [T,sk,~,rd] = BF_ID_Cheby(fun,ri(:),...
                                    ci(:),grid,rk,1e-15,'c',opt);
                            case 'rand'
                                [T,sk,~,rd] = BF_ID_rand(fun,ri(:),...
                                    ci(:),rk,1e-15,'c',tt,opt);
                        end
                    else
                        switch rand_or_cheb
                            case 'cheb'
                                grid = BF_Chey_grid(size(ci(:),1));
                                [T,sk,~,rd] = BF_ID_Cheby(fun,ri(:),...
                                    ci(:),grid,size(ci(:),1),tol,'c',opt);
                            case 'rand'
                                [T,sk,~,rd] = BF_ID_rand(fun,ri(:),...
                                    ci(:),size(ci(:),1),tol,'c',tt,opt);
                        end
                    end
                end
                F.data{l}(b,p).ci  = ci;  % store col block indices
                F.data{l}(b,p).csk = sk;  % store skeleton indices
                F.data{l}(b,p).crd = rd;  % store redundant indices
                T(T==0) = eps;
                F.data{l}(b,p).cT  = T;   % store interpolation matrix
                cidx{b,p} = ci(sk);  % restrict to skeletons
            end
        end
    end
end

% Reform into sparse matrix factors
F.U = cell(L,1);
F.V = cell(L,1);

e = cell((2^Dim),1);
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
for i = 1:(2^Dim)
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
for i = 1:(2^Dim)
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
    np = size(F.data{l},2)/(2^Dim);
    U = cell(np,np);
    V = cell(np,np);
    for q = 1:np  % loop over all partition blocks ...
        for p = 1:np  % ... and compute U and V splitting factors
            [U{p,q},V{p,q}] = BF_IDBF_split_uv_int(F.data{l},p,q,Dim);
        end
    end
    
    % recursively merge U and V factors over all partition blocks
    while np > 1
        np = np/(2^Dim);
        for q = 1:np  % loop over all partition blocks ...
            for p = 1:np  % ... and do 2^Dim * 2^Dim merge until no more
                U{p,q} = BF_IDBF_merge_uv_int(U((2^Dim)*(p-1)+1:...
                    (2^Dim)*p,(2^Dim)*(q-1)+1:(2^Dim)*q),Dim);
                V{p,q} = BF_IDBF_merge_uv_int(V((2^Dim)*(p-1)+1:...
                    (2^Dim)*p,(2^Dim)*(q-1)+1:(2^Dim)*q),Dim);
            end
        end
        U = U(1:np,1:np);
        V = V(1:np,1:np);
    end
    
    % store merged U and V factors
    F.U{l} = U{1};
    F.V{l} = V{1};
end

% form middle sparse S factor
[nb,np] = size(F.data{L});
nbp = nb/np;
S   = cell(np);
tmp = cell(nbp,1);  % work array
for q = 1:np  % loop over secondary partition blocks
    bq = (q-1)*nbp+1:q*nbp;  % blocks in secondary partition
    for p = 1:np  % loop over primary partition blocks
        bp = (p-1)*nbp+1:p*nbp;  % blocks in primary partition
        for b = 1:nbp
            ri = F.data{L}(bp(b),q).ri;
            sk = F.data{L}(bp(b),q).rsk;
            tmp{b} = ri(sk);
        end
        ri = [tmp{:}];  % all remaining rows in (p,q) block
        for b = 1:nbp
            ci = F.data{L}(bq(b),p).ci;
            sk = F.data{L}(bq(b),p).csk;
            tmp{b} = ci(sk);
        end
        ci = [tmp{:}];  % all remaining cols in (p,q) block
        S{p,q} = fun(ri,ci);  % generate skeleton submatrix
    end
end

for l = 1:L
    np = np/(2^Dim);
    for q = 1:np
        for p = 1:np
            S{p,q} = BF_IDBF_merge_s_int(S((2^Dim)*(p-1)+1:(2^Dim)*p,...
                (2^Dim)*(q-1)+1:(2^Dim)*q),Dim);
        end
    end
    S = S(1:np,1:np);
end
F.S = S{1};  % store factorization

% print sparse matrix statistics
nz = nnz(F.S);
for l = 1:L
    nz = nz + nnz(F.U{l}) + nnz(F.V{l});
end

end
