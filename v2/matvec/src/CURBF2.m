function [Factor] = CURBF2(fun_org,xx,kk,nn,r,tol,opt,tag)
% fun_org - function handle representing the BF matrix
% xx - row index set
% kk - column index set
% rk - rank parameter
% tol - set up accuracy
% nn - leaves at least contain nn points
% opt - whether use adaptive rank or fix rank
%       if opt = 0, fix rank rk in the BF no matter what tol is
%       if opt = 1, use a rank less than or equal to rk trying to obtain an
%        accuracy as good as tol; may not be able to achieve tol if rk is
%        too small
% tag - specify how to partition the matrix
% Factor, the butterfly factorization
%
% The CURBF is credited to Eric Michielssen and Amir Boag, MULTILEVEL 
% EVALUATION OF ELECTROMAGNETIC FIELDS FOR THE RAPID SOLUTION OF SCAlTERlNG
% PROBLEMS, Microwave Opt. Technol. Lett., vol. 7, pp. 790-795, Dec. 1994.
%
% Copyright 2018 by Qiyuan Pang, Ken Ho, and Haizhao Yang

if size(xx,1)<size(xx,2), error('xx should be a tall skinny matrix'); end
if size(kk,1)<size(kk,2), error('kk should be a tall skinny matrix'); end
if nargin < 6
    error('more arguments are needed!');
end
if r > 0
    isRk = 1; % use rk in the low-rank approximation
else
    isRk = 0;
end
if tol < 1
    isTol = 1; % use tol in the low-rank approximation
else
    isTol = 0;
end
if nargin <= 7
    opt = 1;
    tag = 'regular';
end
if nargin == 8
    tag = 'regular';
end

[fun,xx,kk,xbox,kbox,npx,npk] = BF_prepbox(fun_org,xx,kk,tag);
% make sure that the smallest dimension of the low-rank submatrices is
% larger than max(8,r)
levels = max(0,min(floor(log2([npx npk])-max(3,ceil(log2(nn))))))+1;
Dim = size(xx,2);

xxset = cell(2*levels+1,1);

xxboxidx = zeros(size(xx));
npxx = npx*2^levels;
xxset{2*levels+1} = cell(prod(npxx),1);
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
    xxset{2*levels+1}{itx} = xx(xxidx(xxIA(itx):xxIA(itx+1)-1),:);
    IIdx = xxidx(xxIA(itx):xxIA(itx+1)-1);
    for i = 1:length(IIdx)
        iidx(ax + i) = IIdx(i);
    end
    ax = ax + length(IIdx);
end

for lvl = 2*levels-1:-1:0
    npxx = npx*2^(lvl-levels);
    xxset{lvl+1} = cell(prod(npxx),1);
    for itx = 1:prod(npxx)
        xxset{lvl+1}{itx} = [];
        for itc = BF_childidx(npxx,itx)
            xxset{lvl+1}{itx} = [ xxset{lvl+1}{itx}; xxset{lvl+2}{itc}];
        end
    end
end

kkset = cell(2*levels+1,1);

kkboxidx = zeros(size(kk));
npkk = npk*2^levels;
kkset{2*levels+1} = cell(prod(npkk),1);
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
    kkset{2*levels+1}{itk} = kk(kkidx(kkIA(itk):kkIA(itk+1)-1),:);
    IIdk = kkidx(kkIA(itk):kkIA(itk+1)-1);
    for i = 1:length(IIdk)
        iidk(ak + i) = IIdk(i);
    end
    ak = ak + length(IIdk);
end

for lvl = 2*levels-1:-1:0
    npkk = npk*2^(lvl-levels);
    kkset{lvl+1} = cell(prod(npkk),1);
    for itk = 1:prod(npkk)
        kkset{lvl+1}{itk} = [];
        for itc = BF_childidx(npkk,itk)
            kkset{lvl+1}{itk} = [ kkset{lvl+1}{itk}; kkset{lvl+2}{itc}];
        end
    end
end

Factor = struct('U',[],'S',[],'V',[]);
Factor.U = cell(1,levels);
Factor.V = cell(1,levels);
N = size(xx,1);

sx = cell(levels,1);
sk = cell(levels,1);

npkk = npk*2^(-levels);
npxx = npx*2^levels;
e = cell(prod(npkk)*prod(npxx),1);
sx{1} = struct('CU',e,'sx1',e,'numr',e);
for i = 1:prod(npkk)
    for j = 1:prod(npxx)
        if isRk & isTol
            [CU,sx1] = CURL(fun,xxset{2*levels+1}{j},kkset{1}{i},r,tol,opt);
        elseif isTol
            [CU,sx1] = CURL(fun,xxset{2*levels+1}{j},kkset{1}{i},size(xxset{2*levels+1}{j},1),tol,opt);
        else
            [CU,sx1] = CURL(fun,xxset{2*levels+1}{j},kkset{1}{i},r,1e-15,opt);
        end
        sx{1}((i-1)*prod(npxx)+j).CU = CU;
        sx{1}((i-1)*prod(npxx)+j).sx1 = sx1;
        sx{1}((i-1)*prod(npxx)+j).numr = size(xxset{2*levels+1}{j},1);
    end
end

kn = prod(npkk);
xn = prod(npxx);
npkk = npk*2^levels;
npxx = npx*2^(-levels);
e = cell(prod(npxx)*prod(npkk),1);
sk{1} = struct('UR',e,'sk1',e,'numc',e);
for k = 1:kn
    for i = 1:prod(npxx)
        x = [];
        for m = 1:xn/prod(npxx)
            x = [x;sx{1}((k-1)*xn+(i-1)*xn/prod(npxx)+m).sx1];
        end
        for j = 1:prod(npkk)/kn
            if isRk & isTol
                [UR,sk1] = CURR(fun,x,kkset{2*levels+1}{(k-1)*prod(npkk)/kn+j},r,tol,opt);
            elseif isTol
                [UR,sk1] = CURR(fun,x,kkset{2*levels+1}{(k-1)*prod(npkk)/kn+j},size(kkset{2*levels+1}{(k-1)*prod(npkk)/kn+j},1),tol,opt);
            else
                [UR,sk1] = CURR(fun,x,kkset{2*levels+1}{(k-1)*prod(npkk)/kn+j},r,1e-15,opt);
            end
            sk{1}((k-1)*prod(npxx)*prod(npkk)/kn+(i-1)*prod(npkk)/kn+j).UR = UR;
            sk{1}((k-1)*prod(npxx)*prod(npkk)/kn+(i-1)*prod(npkk)/kn+j).sk1 = sk1;
            sk{1}((k-1)*prod(npxx)*prod(npkk)/kn+(i-1)*prod(npkk)/kn+j).numc = size(kkset{2*levels+1}{(k-1)*prod(npkk)/kn+j},1);
        end
    end
end


for lvl = 2:levels
    npkk = npk*2^(lvl-levels-1);
    npxx = npx*2^(levels+1-lvl);
    e = cell(prod(npkk)*prod(npxx),1);
    sx{lvl} = struct('CU',e,'sx1',e,'numr',e);
    x1 = prod(npx*2^(levels+2-lvl))/prod(npx*2^(lvl-levels-2));
    k1 = prod(npk*2^(levels+2-lvl))/prod(npk*2^(lvl-levels-2))/2;
    for i = 1:prod(npx*2^(lvl-levels-2))*prod(npk*2^(lvl-levels-2))*2
        k = [];
        for jj = 1:k1
            k = [k;sk{lvl-1}((i-1)*k1+jj).sk1];
        end
        for jj = 1:x1/2
            x = [sx{lvl-1}((((i-mod(i,2))/2+mod(i,2))-1)*x1+(jj-1)*2+1).sx1;sx{lvl-1}((((i-mod(i,2))/2+mod(i,2))-1)*x1+jj*2).sx1];
            if isRk & isTol
                [CU,sx1] = CURL(fun,x,k,r,tol,opt);
            elseif isTol
                [CU,sx1] = CURL(fun,x,k,size(x,1),tol,opt);
            else
                [CU,sx1] = CURL(fun,x,k,r,1e-15,opt);
            end
            sx{lvl}((i-1)*x1/2+jj).CU = CU;
            sx{lvl}((i-1)*x1/2+jj).sx1 = sx1;
            sx{lvl}((i-1)*x1/2+jj).numr = size(x,1);
        end
    end
    
    npkk = npk*2^(levels+1-lvl);
    npxx = npx*2^(lvl-levels-1);
    e = cell(prod(npxx)*prod(npkk),1);
    sk{lvl} = struct('UR',e,'sk1',e,'numc',e);
    x1 = x1/4;
    k1 = k1;
    for i = 1:prod(npx*2^(lvl-levels-2))*prod(npk*2^(lvl-levels-1))*2
        x = [];
        for jj = 1:x1
            x = [x;sx{lvl}((i-1)*x1+jj).sx1];
        end
        for jj = 1:k1/2
            k = [sk{lvl-1}((((i-mod(i,2))/2+mod(i,2))-1)*k1+(jj-1)*2+1).sk1;sk{lvl-1}((((i-mod(i,2))/2+mod(i,2))-1)*k1+jj*2).sk1];
            if isRk & isTol
                [UR,sk1] = CURR(fun,x,k,r,tol,opt);
            elseif isTol
                [UR,sk1] = CURR(fun,x,k,size(k,1),tol,opt);
            else
                [UR,sk1] = CURR(fun,x,k,r,1e-15,opt);
            end
            sk{lvl}((i-1)*k1/2+jj).UR = UR;
            sk{lvl}((i-1)*k1/2+jj).sk1 = sk1;
            sk{lvl}((i-1)*k1/2+jj).numc = size(k,1);
        end
    end
end

npkk = npk*2^(-levels);
npxx = npx*2^levels;
e = cell(prod(npkk),1);
t = struct('XT',e,'YT',e,'ST',e,'totalM',e,'totalN',e);
nnz = 0;
for i = 1:prod(npkk)
    nz = 0;
    for j = 1:prod(npxx)
        c1 = size(sx{1}((i-1)*prod(npxx)+j).sx1,1);
        r1 = sx{1}((i-1)*prod(npxx)+j).numr;
        nz = nz + r1*c1;
    end
    XT = zeros(nz,1);
    YT = zeros(nz,1);
    ST = zeros(nz,1);
    totalM = 0;
    totalN = 0;
    nz = 0;
    for j = 1:prod(npxx)
        CU = sx{1}((i-1)*prod(npxx)+j).CU;
        c1 = size(sx{1}((i-1)*prod(npxx)+j).sx1,1);
        r1 = sx{1}((i-1)*prod(npxx)+j).numr;
        [X Y S] = find(CU);
        ns = length(S);
        XT(nz+(1:ns)) = iidx(totalM + X);
        YT(nz+(1:ns)) = totalN + Y;
        ST(nz+(1:ns)) = S;
        totalM = totalM + r1;
        totalN = totalN + c1;
        nz = nz + ns;
    end
    nnz = nnz + nz;
    t(i).XT = XT(1:nz);
    t(i).YT = YT(1:nz);
    t(i).ST = ST(1:nz);
    t(i).totalM = totalM;
    t(i).totalN = totalN;
end
NN  = 0;
nz = 0;
XT = zeros(nnz,1);
YT = zeros(nnz,1);
ST = zeros(nnz,1);
for i = 1:prod(npkk)
    ns = length(t(i).ST);
    XT(nz+(1:ns)) = t(i).XT;
    YT(nz+(1:ns)) = t(i).YT + NN;
    ST(nz+(1:ns)) = t(i).ST;
    NN  = NN  + t(i).totalN;
    nz = nz + ns;
end
Factor.U{1} = sparse(XT,YT,ST,totalM,NN);

kn = prod(npkk);
xn = prod(npxx);
npkk = npk*2^levels;
npxx = npx*2^(-levels);
e = cell(kn*prod(npxx),1);
t = struct('XT',e,'YT',e,'ST',e,'totalM',e,'totalN',e);
nnz = 0;
fsc = 0;
for k = 1:kn
    for i = 1:prod(npxx)
        nz = 0;
        for j = 1:prod(npkk)/kn
            r1 = size(sk{1}((k-1)*prod(npxx)*prod(npkk)/kn+(i-1)*prod(npkk)/kn+j).sk1,1);
            c1 = sk{1}((k-1)*prod(npxx)*prod(npkk)/kn+(i-1)*prod(npkk)/kn+j).numc;
            nz = nz + c1*r1;
        end
        XT = zeros(nz,1);
        YT = zeros(nz,1);
        ST = zeros(nz,1);
        totalM = 0;
        totalN = 0;
        nz = 0;
        for j = 1:prod(npkk)/kn
            UR = sk{1}((k-1)*prod(npxx)*prod(npkk)/kn+(i-1)*prod(npkk)/kn+j).UR;
            r1 = size(sk{1}((k-1)*prod(npxx)*prod(npkk)/kn+(i-1)*prod(npkk)/kn+j).sk1,1);
            c1 = sk{1}((k-1)*prod(npxx)*prod(npkk)/kn+(i-1)*prod(npkk)/kn+j).numc;
            [X Y S] = find(UR);
            ns = length(S);
            XT(nz+(1:ns)) = totalM + X;
            YT(nz+(1:ns)) = iidk(fsc+totalN + Y);
            ST(nz+(1:ns)) = S;
            totalM = totalM + r1;
            totalN = totalN + c1;
            nz = nz + ns;
        end
        nnz = nnz + nz;
        t((k-1)*prod(npxx)+i).XT = XT(1:nz);
        t((k-1)*prod(npxx)+i).YT = YT(1:nz);
        t((k-1)*prod(npxx)+i).ST = ST(1:nz);
        t((k-1)*prod(npxx)+i).totalM = totalM;
        t((k-1)*prod(npxx)+i).totalN = totalN;
    end
    fsc = fsc + totalN;
end
MM = 0;
nz = 0;
XT = zeros(nnz,1);
YT = zeros(nnz,1);
ST = zeros(nnz,1);
for i = 1:kn*prod(npxx)
    ns = length(t(i).ST);
    XT(nz+(1:ns)) = t(i).XT+MM;
    YT(nz+(1:ns)) = t(i).YT;
    ST(nz+(1:ns)) = t(i).ST;
    MM = MM  + t(i).totalM;
    nz = nz + ns;
end
Factor.V{1} = sparse(XT,YT,ST,MM,fsc);


for lvl = 2:levels
    npkk = npk*2^(lvl-levels-1);
    npxx = npx*2^(levels+1-lvl);
    e = cell(prod(npx*2^(lvl-levels-2))*prod(npk*2^(lvl-levels-2))*2,1);
    t = struct('XT',e,'YT',e,'ST',e,'totalM',e,'totalN',e);
    nnz = 0;
    x1 = prod(npx*2^(levels+2-lvl))/prod(npx*2^(lvl-levels-2));
    k1 = prod(npk*2^(levels+2-lvl))/prod(npk*2^(lvl-levels-2))/2;
    fsr = 0;
    for i = 1:prod(npx*2^(lvl-levels-2))*prod(npk*2^(lvl-levels-2))*2
        nz = 0;
        for jj = 1:x1/2
            r1 = sx{lvl}((i-1)*x1/2+jj).numr;
            c1 = size(sx{lvl}((i-1)*x1/2+jj).sx1,1);
            nz = nz + r1*c1;
        end
        XT = zeros(nz,1);
        YT = zeros(nz,1);
        ST = zeros(nz,1);
        totalM = 0;
        totalN = 0;
        nz = 0;
        for jj = 1:x1/2
            r1 = sx{lvl}((i-1)*x1/2+jj).numr;
            c1 = size(sx{lvl}((i-1)*x1/2+jj).sx1,1);
            CU = sx{lvl}((i-1)*x1/2+jj).CU;
            [X Y S] = find(CU);
            ns = length(S);
            XT(nz+(1:ns)) = fsr + totalM + X;
            YT(nz+(1:ns)) = totalN + Y;
            ST(nz+(1:ns)) = S;
            totalM = totalM + r1;
            totalN = totalN + c1;
            nz = nz + ns;
        end
        nnz = nnz + nz;
        t(i).XT = XT(1:nz);
        t(i).YT = YT(1:nz);
        t(i).ST = ST(1:nz);
        t(i).totalM = totalM;
        t(i).totalN = totalN;
        if  mod(i,2) == 0
            fsr = fsr + totalM;
        end
    end
    NN  = 0;
    nz = 0;
    XT = zeros(nnz,1);
    YT = zeros(nnz,1);
    ST = zeros(nnz,1);
    for i = 1:prod(npx*2^(lvl-levels-2))*prod(npk*2^(lvl-levels-2))*2
        ns = length(t(i).ST);
        XT(nz+(1:ns)) = t(i).XT;
        YT(nz+(1:ns)) = t(i).YT + NN;
        ST(nz+(1:ns)) = t(i).ST;
        NN  = NN  + t(i).totalN;
        nz = nz + ns;
    end
    Factor.U{lvl} = sparse(XT,YT,ST,fsr,NN);
    
    
    npkk = npk*2^(levels+1-lvl);
    npxx = npx*2^(lvl-levels-1);
    e = cell(prod(npx*2^(lvl-levels-2))*prod(npk*2^(lvl-levels-1))*2,1);
    t = struct('XT',e,'YT',e,'ST',e,'totalM',e,'totalN',e);
    x1 = x1/4;
    k1 = k1;
    fsc = 0;
    nnz = 0;
    for i = 1:prod(npx*2^(lvl-levels-2))*prod(npk*2^(lvl-levels-1))*2
        nz = 0;
        for jj = 1:k1/2
            c1 = sk{lvl}((i-1)*k1/2+jj).numc;
            r1 = size(sk{lvl}((i-1)*k1/2+jj).sk1,1);
            nz = nz + c1*r1;
        end
        XT = zeros(nz,1);
        YT = zeros(nz,1);
        ST = zeros(nz,1);
        totalM = 0;
        totalN = 0;
        nz = 0;
        for jj = 1:k1/2
            c1 = sk{lvl}((i-1)*k1/2+jj).numc;
            r1 = size(sk{lvl}((i-1)*k1/2+jj).sk1,1);
            UR = sk{lvl}((i-1)*k1/2+jj).UR;
            [X Y S] = find(UR);
            ns = length(S);
            XT(nz+(1:ns)) = totalM + X;
            YT(nz+(1:ns)) = fsc + totalN +Y;
            ST(nz+(1:ns)) = S;
            totalM = totalM + r1;
            totalN = totalN + c1;
            nz = nz + ns;
        end
        nnz = nnz + nz;
        t(i).XT = XT(1:nz);
        t(i).YT = YT(1:nz);
        t(i).ST = ST(1:nz);
        t(i).totalM = totalM;
        t(i).totalN = totalN;
        if  mod(i,2) == 0
            fsc = fsc + totalN;
        end
    end
    MM = 0;
    nz = 0;
    XT = zeros(nnz,1);
    YT = zeros(nnz,1);
    ST = zeros(nnz,1);
    for i = 1:prod(npx*2^(lvl-levels-2))*prod(npk*2^(lvl-levels-1))*2
        ns = length(t(i).ST);
        XT(nz+(1:ns)) = t(i).XT+MM;
        YT(nz+(1:ns)) = t(i).YT;
        ST(nz+(1:ns)) = t(i).ST;
        MM = MM  + t(i).totalM;
        nz = nz + ns;
    end
    Factor.V{lvl} = sparse(XT,YT,ST,MM,fsc);
end


x1 = prod(npx*2)/prod(npx*2^(-1));
k1 = prod(npk*2)/prod(npk*2^(-1));
e = cell(prod(npx)*prod(npk)/4,1);
ss = struct('x',e,'k',e,'r',e,'c',e);
nz = 0;
for i = 1:prod(npx)*prod(npk)/4
    x = [];
    for j = 1:x1
        x = [x;sx{levels}((i-1)*x1+j).sx1];
    end
    k = [];
    for j = 1:k1
        k = [k;sk{levels}((i-1)*k1+j).sk1];
    end
    ss(i).x = x;
    ss(i).k = k;
    ss(i).r = size(x,1);
    ss(i).c = size(k,1);
    nz = nz + ss(i).r*ss(i).c;
end
XT = zeros(nz,1);
YT = zeros(nz,1);
ST = zeros(nz,1);
nz = 0;
totalM = 0;
totalN = 0;
for i = 1:prod(npx)*prod(npk)/4
    C = fun(ss(i).x,ss(i).k);
    [X Y S] = find(C);
    ns = length(S);
    XT(nz+(1:ns)) = totalM + X;
    YT(nz+(1:ns)) = totalN + Y;
    ST(nz+(1:ns)) = S;
    totalM = totalM + ss(i).r;
    totalN = totalN + ss(i).c;
    nz = nz + ns;
end
Factor.S = sparse(XT,YT,ST,totalM,totalN);
end