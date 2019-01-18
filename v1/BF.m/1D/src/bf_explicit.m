function Factor = bf_explicit(fun, xx, xbox, kk, kbox, mR, tol, disp_flag)

if(disp_flag)
    fprintf('Butterfly Factorization explicit version started...\n\n');
end

% Nx is the square root of the number of target points in space
[Nx,~] = size(xx);
% Np is the square root of the number of source points in phase
[Nk,~] = size(kk);

tR=3*mR;

% npx is the number of blocks of each dimension in space
npx = 2^ceil(log2(sqrt(Nx))+0.5);
% npk is the number of blocks of each dimension in phase
npk = 2^ceil(log2(sqrt(Nk))+0.5);

if(disp_flag)
    fprintf('Space  Frequence\n');
    fprintf('  npx    npk\n');
    fprintf(' %4d   %4d\n\n',npx,npk);
end

P = cell(npx,npk);
Ridx = cell(npx,npk);
Cidx = cell(npx,npk);
rs = cell(npx,npk);
cs = cell(npx,npk);

xidx = bf_prep(xx,xbox,npx);
kidx = bf_prep(kk,kbox,npk);

levels = max(ceil(log2(Nx/npx/mR/2)),0);
LS = 2*mR^2*npk*npx;

if(disp_flag)
    fprintf('Compression levels: %d\n',levels);
    fprintf('Preallocated sparse matrix size: %d, about %.2f GB\n', ...
        LS,LS*(levels+2)*(2*8+16)/1024/1024/1024);
end

CPreSpr = repmat(struct('XT',zeros(LS,1),'YT',zeros(LS,1), ...
    'ST',zeros(LS,1),'Height',0,'Width',0,'Offset',0),levels,1);
WPreSpr = struct('XT',zeros(2*LS,1),'YT',zeros(2*LS,1), ...
    'ST',zeros(2*LS,1),'Height',Nx,'Width',0,'Offset',0);

for x=1:npx
    U = cell(npk,1);
    for k=1:npk
        ik = kidx{k};
        ix = xidx{x};
        [U{k},Stmp,~,Ridx{x,k},Cidx{x,k},rs{x,k},cs{x,k}] = lowrank(xx(ix), kk(ik), fun, tol, tR, mR);
        P{x,k} = 1./diag(Stmp);
        U{k} = U{k}*Stmp;
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
    [WPreSpr,CPreSpr] = compression1D(WPreSpr,CPreSpr,U,xxsub,xidx{x},xsubbox,mR,tol,1,levels);
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
ATol = cell(levels,1);
for l=1:levels
    Aid = find(CPreSpr(l).ST~=0);
    ATol{l} = sparse(CPreSpr(l).XT(Aid),CPreSpr(l).YT(Aid),CPreSpr(l).ST(Aid));
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
    memsize = whos('ATol');
    fprintf('Compressed A Memory: %.2f GB\n',memsize.bytes/1024/1024/1024);
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
        [~,Stmp,V{x}] = lowrankidx(xx(ix),kk(ik),fun,tol,mR,Ridx{x,k},Cidx{x,k},rs{x,k},cs{x,k});
        V{x} = V{x}*Stmp;
    end
    kksub = kk(kidx{k});
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
    [WPreSpr,CPreSpr] = compression1D(WPreSpr,CPreSpr,V,kksub,kidx{k},ksubbox,mR,tol,1,levels);
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
BTol = cell(levels,1);
for l=1:levels
    Bid = find(CPreSpr(l).ST~=0);
    BTol{l} = sparse(CPreSpr(l).XT(Bid),CPreSpr(l).YT(Bid),CPreSpr(l).ST(Bid));
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
    memsize = whos('BTol');
    fprintf('Compressed B Memory: %.2f GB\n',memsize.bytes/1024/1024/1024);
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
SigmaM = sparse(XT,YT,ST);
clear XT YT ST;
if(disp_flag)
    memsize = whos('SigmaM');
    fprintf('Compressed M Memory: %.2f GB\n',memsize.bytes/1024/1024/1024);
    clear memsize;
end

Factor = struct('U',[],'ATol',[],'SigmaM',[],'BTol',[],'V',[]);
Factor.U = USpr;
clear USpr;
Factor.ATol = ATol;
clear ATol;
Factor.SigmaM = SigmaM;
clear SigmaM;
Factor.BTol = BTol;
clear BTol;
Factor.V = VSpr;
clear VSpr;

end
