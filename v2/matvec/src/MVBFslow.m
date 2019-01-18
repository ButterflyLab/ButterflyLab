function Factor = MVBFslow(fun, fun_adj, xx, xbox, kk, kbox, NG, tol, disp_flag)
% O(N^1.5 log N) operation complexity and O(N^1.5) memory complexity.
%
% Reference: Y. Li, H. Yang, E. Martin, K. Ho and L. Ying, Butterfly 
% Factorization, SIAM Multiscale Modeling and Simulation, 2015.


if size(xx,1)<size(xx,2), error('xx should be a tall skinny matrix'); end
if size(kk,1)<size(kk,2), error('kk should be a tall skinny matrix'); end

if(disp_flag)
    fprintf('Butterfly Factorization implicit version started...\n\n');
end

% Nx is the square root of the number of target points in space
[Nx,~] = size(xx);
% Nk is the square root of the number of source points in phase
[Nk,~] = size(kk);

tR=NG+5;

% npx is the number of blocks of each dimension in space
npx = 2^ceil(log2(sqrt(Nx))+0.5);
% npk is the number of blocks of each dimension in phase
npk = 2^ceil(log2(sqrt(Nk))+0.5);

if(disp_flag)
    fprintf('Space  Frequence\n');
    fprintf('  npx    npk\n');
    fprintf(' %4d   %4d\n\n',npx,npk);
end

U_uc = cell(npk,npx);
V_uc = cell(npk,npx);
P = cell(npk,npx);

xidx = BF_prep(xx,xbox,npx);
kidx = BF_prep(kk,kbox,npk);

f1all = randn(Nk,tR);% + sqrt(-1)*randn(Nk,tR);
f2all = randn(Nx,tR);% + sqrt(-1)*randn(Nx,tR);


% make sure that the smallest dimension of the low-rank submatrices is
% larger than max(8,NG)
levels = max(0,min(floor(log2([npx npk])-max(3,ceil(log2(NG))))));  

if(disp_flag)
    fprintf('Compression levels: %d\n',levels);
    fprintf('Estimated Uncompressed U,V matrix size: %d, about %.2f GB\n\n', ...
        Nx*npk*tR,2*Nx*npk*tR*16/1024/1024/1024);
end

if(disp_flag)
    t_start = cputime;
end

for k=1:npk
    ik = kidx{k};
    f1 = zeros(Nk,tR);
    f1(ik,:)=f1all(ik,:);
    BR = fun(f1);
    for x=1:npx
        ix = xidx{x};
        U_uc{x,k}=BR(ix,:);
    end
    if(disp_flag)
        if( k==1 )
            fprintf('Multiply U block: (%4d/%4d)',k,npk);
        else
            fprintf('\b\b\b\b\b\b\b\b\b\b\b');
            fprintf('(%4d/%4d)',k,npk);
        end
    end
end

if(disp_flag)
    fprintf('\n');
    fprintf('Multiply for U block cost %.2f seconds\n',cputime - t_start);
    memsize = whos('U_uc');
    fprintf('Uncompressed U Memory: %.2f GB\n\n',memsize.bytes/1024/1024/1024);
    clear memsize;
end
clear BR f1;

if(disp_flag)
    t_start = cputime;
end

for x=1:npx
    ix = xidx{x};
    f2 = zeros(Nx,tR);
    f2(ix,:)=f2all(ix,:);
    BHR = fun_adj(f2);
    for k=1:npk
        ik = kidx{k};
        V_uc{x,k}=BHR(ik,:);
    end
    if(disp_flag)
        if( x==1 )
            fprintf('Multiply V block: (%4d/%4d)',x,npx);
        else
            fprintf('\b\b\b\b\b\b\b\b\b\b\b');
            fprintf('(%4d/%4d)',x,npx);
        end
    end
end

if(disp_flag)
    fprintf('\n');
    fprintf('Multiply for V block cost %.2f seconds\n',cputime - t_start);
    memsize = whos('V_uc');
    fprintf('Uncompressed V Memory: %.2f GB\n\n',memsize.bytes/1024/1024/1024);
    clear memsize;
end
clear BHR f2;

LS = 2*NG^2*npk*npx;
CPreSpr = repmat(struct('XT',zeros(LS,1),'YT',zeros(LS,1), ...
    'ST',zeros(LS,1),'Height',0,'Width',0,'Offset',0),levels,1);
WPreSpr = struct('XT',zeros(2*LS,1),'YT',zeros(2*LS,1), ...
    'ST',zeros(2*LS,1),'Height',Nx^2,'Width',0,'Offset',0);

for x=1:npx
    U = cell(npk,1);
    for k=1:npk
        ik = kidx{k};
        ix = xidx{x};
        f1=f1all(ik,:);
        f2=f2all(ix,:);
        BR = U_uc{x,k};
        BHR = V_uc{x,k};

        [VC,~] = qr(BR,0);
        [VR,~] = qr(BHR,0);

        RrVC = f2'*VC;
        VRRc = VR'*f1;
        [Utmp,Stmp,~] = BF_svdtrunc_rank(pinv(RrVC) * f2'*BR * pinv(VRRc),NG,tol);
        U{k} = VC*Utmp*Stmp;
        P{x,k} = 1./diag(Stmp);
    end
    xxsub = xx(xidx{x},:);
    xos = xbox(1);
    xlen = xbox(2)-xbox(1);
    xsubbox = [ xos+(x-1)*xlen/npx, xos+x*xlen/npx ];
    if(disp_flag)
        if( x==1 )
            fprintf('Compress U block: (%4d/%4d)',x,npx);
        else
            fprintf('\b\b\b\b\b\b\b\b\b\b\b');
            fprintf('(%4d/%4d)',x,npx);
        end
    end
    [WPreSpr,CPreSpr] = BF_compression1D(WPreSpr,CPreSpr,U,xxsub,xidx{x},xsubbox,NG,tol,1,levels);
end

if(disp_flag)
    fprintf('\n');
end

Uid = find(WPreSpr.ST~=0);
USpr = sparse(WPreSpr.XT(Uid),WPreSpr.YT(Uid),WPreSpr.ST(Uid));
if(disp_flag)
    if(length(WPreSpr.XT)>2*LS)
        fprintf('Bad preallocation U, %d is required\n',length(WPreSpr.XT));
    end
end
clear Uid;
GTol = cell(levels,1);
for l=1:levels
    Aid = find(CPreSpr(l).ST~=0);
    GTol{l} = sparse(CPreSpr(l).XT(Aid),CPreSpr(l).YT(Aid),CPreSpr(l).ST(Aid));
    if(disp_flag)
        if(length(CPreSpr(l).XT)>LS)
            fprintf('Bad preallocation A, %d is required\n',length(CPreSpr(l).XT));
        end
    end
end
clear Aid;

if(disp_flag)
    memsize = whos('USpr');
    fprintf('Compressed U Memory: %.2f GB\n',memsize.bytes/1024/1024/1024);
    memsize = whos('GTol');
    fprintf('Compressed A Memory: %.2f GB\n\n',memsize.bytes/1024/1024/1024);
    clear memsize;
end

CPreSpr = repmat(struct('XT',zeros(LS,1),'YT',zeros(LS,1), ...
    'ST',zeros(LS,1),'Height',0,'Width',0,'Offset',0),levels,1);
WPreSpr = struct('XT',zeros(2*LS,1),'YT',zeros(2*LS,1), ...
    'ST',zeros(2*LS,1),'Height',Nk,'Width',0,'Offset',0);

for k=1:npk
    V = cell(npx,1);
    for x=1:npx
        ik = kidx{k};
        ix = xidx{x};
        f1=f1all(ik,:);
        f2=f2all(ix,:);
        BR = U_uc{x,k};
        BHR = V_uc{x,k};

        [VC,~] = qr(BR,0);
        [VR,~] = qr(BHR,0);

        RrVC = f2'*VC;
        VRRc = VR'*f1;
        [~,Stmp,Vtmp] = BF_svdtrunc_rank(pinv(RrVC) * f2'*BR * pinv(VRRc),NG,tol);
        V{x} = VR*Vtmp*Stmp;
    end
    kksub = kk(kidx{k},:);
    kos = kbox(1);
    klen = kbox(2)-kbox(1);
    ksubbox = [ kos+(k-1)*klen/npk, kos+k*klen/npk ];
    if(disp_flag)
        if( k==1 )
            fprintf('Compress V block: (%4d/%4d)',k,npk);
        else
            fprintf('\b\b\b\b\b\b\b\b\b\b\b');
            fprintf('(%4d/%4d)',k,npk);
        end
    end
    [WPreSpr,CPreSpr] = BF_compression1D(WPreSpr,CPreSpr,V,kksub,kidx{k},ksubbox,NG,tol,1,levels);
end

if(disp_flag)
    fprintf('\n');
end

Vid = find(WPreSpr.ST~=0);
VSpr = sparse(WPreSpr.XT(Vid),WPreSpr.YT(Vid),WPreSpr.ST(Vid));
if(disp_flag)
    if(length(WPreSpr.XT)>2*LS)
        fprintf('Bad preallocation V, %d is required\n',length(WPreSpr.XT));
    end
end
clear Vid;
HTol = cell(levels,1);
for l=1:levels
    Bid = find(CPreSpr(l).ST~=0);
    HTol{l} = sparse(CPreSpr(l).XT(Bid),CPreSpr(l).YT(Bid),CPreSpr(l).ST(Bid));
    if(disp_flag)
        if(length(CPreSpr(l).XT)>LS)
            fprintf('Bad preallocation B, %d is required\n',length(CPreSpr(l).XT));
        end
    end
end
clear Bid;
clear WPreSpr CPreSpr;
if(disp_flag)
    memsize = whos('VSpr');
    fprintf('Compressed V Memory: %.2f GB\n',memsize.bytes/1024/1024/1024);
    memsize = whos('HTol');
    fprintf('Compressed B Memory: %.2f GB\n\n',memsize.bytes/1024/1024/1024);
    clear memsize;
end

totalH = zeros(npx,1);
for x=1:npx
    for k=1:npk
        totalH(x) = totalH(x) + length(P{x,k});
    end
end
currentH = zeros(npx,1);
for x=1:npx
    if(x>1)
        currentH(x) = currentH(x-1)+totalH(x-1);
    end
end

totalW = zeros(npk,1);
for k=1:npk
    for x=1:npx
        totalW(k) = totalW(k) + length(P{x,k});
    end
end
currentW = zeros(npk,1);
for k=1:npk
    if(k>1)
        currentW(k) = currentW(k-1)+totalW(k-1);
    end
end

totalel = 0;
for x=1:npx
    for k=1:npk
        totalel = totalel + length(P{x,k});
    end
end

offset = 0;
XT = zeros(totalel,1);
YT = zeros(totalel,1);
ST = zeros(totalel,1);
for x=1:npx
    localH = currentH(x);
    for k=1:npk
        Mlen = length(P{x,k});
        X = localH+(1:Mlen);
        Y = currentW(k)+(1:Mlen);
        idx = offset+1:offset+numel(X);
        XT(idx) = X(:);
        YT(idx) = Y(:);
        ST(idx) = P{x,k};
        if(~isempty(idx))
            offset = idx(end);
        end
        localH = localH + Mlen;
        currentW(k) = currentW(k) + Mlen;
    end
end
M = sparse(XT,YT,ST);
clear XT YT ST;
if(disp_flag)
    memsize = whos('M');
    fprintf('Compressed M Memory: %.2f GB\n\n',memsize.bytes/1024/1024/1024);
    clear memsize;
end

Factor = struct('U',[],'GTol',[],'M',[],'HTol',[],'V',[]);
Factor.U = USpr;
clear USpr;
Factor.GTol = GTol;
clear GTol;
Factor.M = M;
clear M;
Factor.HTol = HTol;
clear HTol;
Factor.V = VSpr';
clear VSpr;

Factor = BF_resize(Factor, Nx,Nk);

end
