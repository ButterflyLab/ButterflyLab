function [Factor,timeinfo] = MVRSslow_polar2D(N, fun, fun_adj, xx, xbox, kk, kbox, mR, tol, disp_flag, runlimit)
% O(N^1.5 log N) operation complexity and O(N^1.5) memory complexity.
%
% Reference: Y. Li, H. Yang, E. Martin, K. Ho and L. Ying, Butterfly 
% Factorization, SIAM Multiscale Modeling and Simulation, 2015.


if size(xx,1)<size(xx,2), error('xx should be a tall skinny matrix'); end
if size(kk,1)<size(kk,2), error('kk should be a tall skinny matrix'); end

if(nargin<11)
    runlimit = -1;
end

if(runlimit>0)
    savefile = ['tmp/pbf_implicit_' num2str(N) '_' num2str(mR) '.mat'];
    if(~exist('tmp/', 'dir'))
        mkdir('tmp/');
    end
end

if(disp_flag)
    fprintf('Butterfly Factorization implicit version started...\n\n');
end

pp = pbf_k2p(N,kk);
pbox = (kbox+N/2)/N;
pbox(1,2) = pbox(1,2)+1/N;

% Nx is the square root of the number of target points in space
[Nxx,~] = size(xx);
Nx = floor(sqrt(Nxx));
% Np is the square root of the number of source points in phase
[Npp,~] = size(pp);
Np = floor(sqrt(Npp));

tR=mR+5;

% npx is the number of blocks of each dimension in space
npx1 = 2^ceil(log2(sqrt(Nx)));
npx2 = 2^ceil(log2(sqrt(Nx)));
% npp is the number of blocks of each dimension in phase
npp1 = 2^ceil(log2(sqrt(Np)));
npp2 = 4*2^ceil(log2(sqrt(Np)));

if(disp_flag)
    fprintf('   Space      Frequence\n');
    fprintf('  npx1  npx2  npp1  npp2\n');
    fprintf('  %4d  %4d  %4d  %4d\n\n',npx1,npx2,npp1,npp2);
end

xidx = bf_prep(xx,xbox,npx1,npx2);
pidx = bf_prep(pp,pbox,npp1,npp2);

if(runlimit>0 && exist(savefile,'file'))
    tmpdata = load(savefile);
    p1s = tmpdata.p1s;
    x1s = tmpdata.x1s;
    f1all = tmpdata.f1all;
    f2all = tmpdata.f2all;
    U_uc = tmpdata.U_uc;
    V_uc = tmpdata.V_uc;
    all_t = tmpdata.all_t;
    clear tmpdata;
else
    p1s = 1;
    x1s = 1;
    f1all = randn(Npp,tR);% + 1i*randn(Nkk,tR);
    f2all = randn(Nxx,tR);% + 1i*randn(Nxx,tR);
    U_uc = cell(npx1,npx2,npp1,npp2);
    V_uc = cell(npx1,npx2,npp1,npp2);
    all_t = 0;
end

all_t_start = cputime;

XYmesh = cell(mR,mR,2);
for Xiter = 1:mR
    for Yiter = 1:mR
        [X,Y] = meshgrid(1:Xiter,1:Yiter);
        XYmesh{Xiter,Yiter,1} = X';
        XYmesh{Xiter,Yiter,2} = Y';
    end
end

levels = max(ceil(log2(Nxx/npx1/npx2/mR/4)),0);
LS = 4*mR^2*npp1*npp2*npx1*npx2;

if(disp_flag)
    fprintf('Compression levels: %d\n',levels);
    fprintf('Preallocated sparse matrix size: %d, about %.2f GB\n\n', ...
        LS,LS*(levels+2)*(2*8+16)/1024/1024/1024);
end

P = cell(npx1,npx2,npp1,npp2);
Vtmpcell = cell(npx1,npx2,npp1,npp2);

if(disp_flag)
    t_start = cputime;
end

for p1=p1s:npp1
    for p2=1:npp2
        ip = pidx{p1,p2};
        f1 = zeros(Npp,tR);
        f1(ip,:)=f1all(ip,:);
        BR = fun(f1);
        for x1=1:npx1
            for x2=1:npx2
                ix = xidx{x1,x2};
                U_uc{x1,x2,p1,p2}=BR(ix,:);
            end
        end
        if(disp_flag)
            if( p1==p1s && p2==1 )
                fprintf('Multiply U block: (%4d/%4d) by (%4d/%4d)',p1,npp1,p2,npp2);
            else
                fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b');
                fprintf('(%4d/%4d) by (%4d/%4d)',p1,npp1,p2,npp2);
            end
        end
    end
    if( runlimit>0 && cputime-all_t_start>runlimit*60*60 )
        fprintf('\n');
        p1s = p1+1;
        x1s = 1;
        all_t = all_t + cputime - all_t_start;
        save(savefile,'p1s','x1s','f1all','f2all','U_uc','V_uc','all_t','-v7.3');
        timeinfo = -1;
        Factor = [];
        return;
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

for x1=x1s:npx1
    for x2=1:npx2
        ix = xidx{x1,x2};
        f2 = zeros(Nxx,tR);
        f2(ix,:)=f2all(ix,:);
        BHR = fun_adj(f2);
        for p1=1:npp1
            for p2=1:npp2
                ip = pidx{p1,p2};
                V_uc{x1,x2,p1,p2}=BHR(ip,:);
            end
        end
        if(disp_flag)
            if( x1==x1s && x2==1 )
                fprintf('Multiply V block: (%4d/%4d) by (%4d/%4d)',x1,npx1,x2,npx2);
            else
                fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b');
                fprintf('(%4d/%4d) by (%4d/%4d)',x1,npx1,x2,npx2);
            end
        end
    end
    if( runlimit>0 && cputime-all_t_start > runlimit*60*60 )
        fprintf('\n');
        p1s = npp1+1;
        x1s = x1+1;
        all_t = all_t + cputime - all_t_start;
        save(savefile,'p1s','x1s','f1all','f2all','U_uc','V_uc','all_t','-v7.3');
        timeinfo = -1;
        Factor = [];
        return;
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

CPreSpr = repmat(struct('XT',zeros(LS,1),'YT',zeros(LS,1), ...
    'ST',zeros(LS,1),'Height',0,'Width',0,'Offset',0),levels,1);
WPreSpr = struct('XT',zeros(2*LS,1),'YT',zeros(2*LS,1), ...
    'ST',zeros(2*LS,1),'Height',Nx^2,'Width',0,'Offset',0);

for x1=1:npx1
    for x2=1:npx2
        U = cell(npp1,npp2);
        for p1=1:npp1
            for p2=1:npp2
                ip = pidx{p1,p2};
                ix = xidx{x1,x2};
                f1=f1all(ip,:);
                f2=f2all(ix,:);
                BR = U_uc{x1,x2,p1,p2};
                BHR = V_uc{x1,x2,p1,p2};

                [VC,~] = qr(BR,0);
                [VR,~] = qr(BHR,0);

                RrVC = f2'*VC;
                VRRc = VR'*f1;
                [Utmp,Stmp,Vtmp] = svdrand(pinv(RrVC) * f2'*BR * pinv(VRRc),mR,tol);
                U{p1,p2} = VC*Utmp*Stmp;
                P{x1,x2,p1,p2} = 1./diag(Stmp);
		Vtmpcell{x1,x2,p1,p2} = Vtmp*Stmp;
            end
        end
        xxsub = xx(xidx{x1,x2},:);
        x1os = xbox(1,1);
        x1len = xbox(1,2)-xbox(1,1);
        x2os = xbox(2,1);
        x2len = xbox(2,2)-xbox(2,1);
        xsubbox = [ x1os+(x1-1)*x1len/npx1, x1os+x1*x1len/npx1; ...
            x2os+(x2-1)*x2len/npx2, x2os+x2*x2len/npx2];
        if(disp_flag)
            if( x1==1 && x2==1 )
                fprintf('Compress U block: (%4d/%4d) by (%4d/%4d)',x1,npx1,x2,npx2);
            else
                fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b');
                fprintf('(%4d/%4d) by (%4d/%4d)',x1,npx1,x2,npx2);
            end
        end
        [WPreSpr,CPreSpr] = compression2D(WPreSpr,CPreSpr,U,xxsub,xidx{x1,x2},xsubbox,mR,tol,1,levels,XYmesh);
    end
end

clear U_uc;

if(disp_flag)
    fprintf('\n');
end

Uid = find(WPreSpr.ST~=0);
USpr = sparse(WPreSpr.XT(Uid),WPreSpr.YT(Uid),WPreSpr.ST(Uid));
clear Uid;
ATol = cell(levels,1);
for l=1:levels
    Aid = find(CPreSpr(l).ST~=0);
    ATol{l} = sparse(CPreSpr(l).XT(Aid),CPreSpr(l).YT(Aid),CPreSpr(l).ST(Aid));
end
clear Aid;

if(disp_flag)
    memsize = whos('USpr');
    fprintf('Compressed U Memory: %.2f GB\n',memsize.bytes/1024/1024/1024);
    memsize = whos('ATol');
    fprintf('Compressed A Memory: %.2f GB\n\n',memsize.bytes/1024/1024/1024);
    clear memsize;
end

CPreSpr = repmat(struct('XT',zeros(LS,1),'YT',zeros(LS,1), ...
    'ST',zeros(LS,1),'Height',0,'Width',0,'Offset',0),levels,1);
WPreSpr = struct('XT',zeros(2*LS,1),'YT',zeros(2*LS,1), ...
    'ST',zeros(2*LS,1),'Height',Npp,'Width',0,'Offset',0);

for p1=1:npp1
    for p2=1:npp2
        V = cell(npx1,npx2);
        for x1=1:npx1
            for x2=1:npx2
                BHR = V_uc{x1,x2,p1,p2};

                [VR,~] = qr(BHR,0);

                V{x1,x2} = VR*Vtmpcell{x1,x2,p1,p2};
            end
        end
        ppsub = pp(pidx{p1,p2},:);
        p1os = pbox(1,1);
        p1len = pbox(1,2)-pbox(1,1);
        p2os = pbox(2,1);
        p2len = pbox(2,2)-pbox(2,1);
        psubbox = [ p1os+(p1-1)*p1len/npp1, p1os+p1*p1len/npp1; ...
            p2os+(p2-1)*p2len/npp2, p2os+p2*p2len/npp2];
        if(disp_flag)
            if( p1==1 && p2==1 )
                fprintf('Compress V block: (%4d/%4d) by (%4d/%4d)',p1,npp1,p2,npp2);
            else
                fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b');
                fprintf('(%4d/%4d) by (%4d/%4d)',p1,npp1,p2,npp2);
            end
        end
        [WPreSpr,CPreSpr] = compression2D(WPreSpr,CPreSpr,V,ppsub,pidx{p1,p2},psubbox,mR,tol,1,levels,XYmesh);
    end
end

clear V_uc Vtmpcell;

if(disp_flag)
    fprintf('\n');
end

Vid = find(WPreSpr.ST~=0);
VSpr = sparse(WPreSpr.XT(Vid),WPreSpr.YT(Vid),WPreSpr.ST(Vid));
clear Vid;
BTol = cell(levels,1);
for l=1:levels
    Bid = find(CPreSpr(l).ST~=0);
    BTol{l} = sparse(CPreSpr(l).XT(Bid),CPreSpr(l).YT(Bid),CPreSpr(l).ST(Bid));
end
clear Bid;
clear WPreSpr CPreSpr;
if(disp_flag)
    memsize = whos('VSpr');
    fprintf('Compressed V Memory: %.2f GB\n',memsize.bytes/1024/1024/1024);
    memsize = whos('BTol');
    fprintf('Compressed B Memory: %.2f GB\n\n',memsize.bytes/1024/1024/1024);
    clear memsize;
end

totalH = zeros(npx1*npx2,1);
for x1=1:npx1
    for x2=1:npx2
        x = (x1-1)*npx2+x2;
        for p1=1:npp1
            for p2=1:npp2
                totalH(x) = totalH(x) + length(P{x1,x2,p1,p2});
            end
        end
    end
end
currentH = zeros(npx1*npx2,1);
for x1=1:npx1
    for x2 =1:npx2
        x = (x1-1)*npx2+x2;
        if(x>1)
            currentH(x) = currentH(x-1)+totalH(x-1);
        end
    end
end

totalW = zeros(npp1*npp2,1);
for p1=1:npp1
    for p2=1:npp2
        p = (p1-1)*npp2+p2;
        for x1=1:npx1
            for x2=1:npx2
                totalW(p) = totalW(p) + length(P{x1,x2,p1,p2});
            end
        end
    end
end
currentW = zeros(npp1*npp2,1);
for p1=1:npp1
    for p2=1:npp2
        p = (p1-1)*npp2+p2;
        if(p>1)
            currentW(p) = currentW(p-1)+totalW(p-1);
        end
    end
end

totalel = 0;
for x1=1:npx1
    for x2=1:npx2
        for p1=1:npp1
            for p2=1:npp2
                totalel = totalel + length(P{x1,x2,p1,p2});
            end
        end
    end
end

offset = 0;
XT = zeros(totalel,1);
YT = zeros(totalel,1);
ST = zeros(totalel,1);
for x1=1:npx1
    for x2=1:npx2
        x = (x1-1)*npx2+x2;
        localH = currentH(x);
        for p1=1:npp1
            for p2=1:npp2
                p = (p1-1)*npp2+p2;
                Mlen = length(P{x1,x2,p1,p2});
                X = localH+(1:Mlen);
                Y = currentW(p)+(1:Mlen);
                idx = offset+1:offset+numel(X);
                XT(idx) = X(:);
                YT(idx) = Y(:);
                ST(idx) = P{x1,x2,p1,p2};
                if(~isempty(idx))
                    offset = idx(end);
                end
                localH = localH + Mlen;
                currentW(p) = currentW(p) + Mlen;
            end
        end
    end
end
SigmaM = sparse(XT,YT,ST);
clear XT YT ST;
if(disp_flag)
    memsize = whos('SigmaM');
    fprintf('Compressed M Memory: %.2f GB\n\n',memsize.bytes/1024/1024/1024);
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
timeinfo = all_t + cputime - all_t_start;

end
