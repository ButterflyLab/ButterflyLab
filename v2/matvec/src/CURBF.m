function [Factor,Rcomp] = CURBF(fun_org,xx,kk,r,tol,tag)
% This is the main file to compute the CUR butterfly factorization.
% The CURBF is credited to Eric Michielssen and Amir Boag, MULTILEVEL 
% EVALUATION OF ELECTROMAGNETIC FIELDS FOR THE RAPID SOLUTION OF SCAlTERlNG
% PROBLEMS, Microwave Opt. Technol. Lett., vol. 7, pp. 790-795, Dec. 1994.
%
% Copyright 2017 by Yingzhou Li and Haizhao Yang

% xx = xx(:); kk = kk(:);
if size(xx,1)<size(xx,2), error('xx should be a tall skinny matrix'); end
if size(kk,1)<size(kk,2), error('kk should be a tall skinny matrix'); end
if nargin == 5
    tag = 'Regular';
end

[fun,xx,kk,xbox,kbox,npx,npk] = BF_prepbox(fun_org,xx,kk,tag);

% make sure that the smallest dimension of the low-rank submatrices is
% larger than max(8,r)
levels = max(0,min(floor(log2([npx npk])-max(3,ceil(log2(r)))))); 
Dim = size(xx,2);

xxset = cell(levels+1,1);

xxboxidx = zeros(size(xx));
npxx = npx*2^levels;
xxset{levels+1} = cell(prod(npxx),1);
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
for itx = 1:prod(npxx)
    xxset{levels+1}{itx} = xx(xxidx(xxIA(itx):xxIA(itx+1)-1),:);
end

for lvl = levels-1:-1:0
    npxx = npx*2^lvl;
    xxset{lvl+1} = cell(prod(npxx),1);
    for itx = 1:prod(npxx)
        xxset{lvl+1}{itx} = [];
        for itc = BF_childidx(npxx,itx)
            xxset{lvl+1}{itx} = [ xxset{lvl+1}{itx}; xxset{lvl+2}{itc}];
        end
    end
end

kkset = cell(levels+1,1);

kkboxidx = zeros(size(kk));
npkk = npk*2^levels;
kkset{levels+1} = cell(prod(npkk),1);
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
for itk = 1:prod(npkk)
    kkset{levels+1}{itk} = kk(kkidx(kkIA(itk):kkIA(itk+1)-1),:);
end

for lvl = levels-1:-1:0
    npkk = npk*2^lvl;
    kkset{lvl+1} = cell(prod(npkk),1);
    for itk = 1:prod(npkk)
        kkset{lvl+1}{itk} = [];
        for itc = BF_childidx(npkk,itk)
            kkset{lvl+1}{itk} = [ kkset{lvl+1}{itk}; kkset{lvl+2}{itc}];
        end
    end
end

%===============================================================
% The butterfly factorization is of the form
%  K = U Gl ... G2 G1 M H1 H2 ... Hl V
% where U is stored in 'U', G1 G2 ... Gl are stored in 'GTol',
% M is stored in 'M', H1 H2 ... Hl are stored in 'HTol',
% and V is stored in 'V'.
%===============================================================
Factor = struct('U',[],'GTol',[],'M',[],'HTol',[],'V',[]);

%---------------------------------------------------------------
%   M construction

Mcell = cell(prod(npx),prod(npk));
sx = cell(prod(npk),prod(npx));
sk = cell(prod(npx),prod(npk));

for itx = 1:prod(npx)
    xxval = xxset{1}{itx};
    for itk = 1:prod(npk)
        kkval = kkset{1}{itk};
        [ Mcell{itx,itk}, sx{itk,itx}, sk{itx,itk} ] = CUR(fun, ...
            xxval, kkval, r, tol);
        
        
    end
end

%---------------------------------------------------------------
%   G construction

GTolcell = cell(levels,1);
for ell = 1:levels
    npxx = npx*2^ell;
    npkk = npk/2^ell;
    skpre = sk;
    sk = cell(prod(npxx),prod(npkk));
    
    GTolcell{ell} = cell(prod(npxx),prod(npkk)*2^Dim);
    for itx = 1:prod(npxx)
        xxval = xxset{ell+1}{itx};
        xxparidx = BF_parentidx(npxx,itx);
        for itk = 1:prod(npkk)
            cursk = [];
            for itk_child = BF_childidx(npkk,itk)
                cursk = [cursk; skpre{xxparidx,itk_child}];
            end
            [ UR, sk{itx,itk} ] = CURR(fun, xxval, cursk, r, tol);
            offset = 0;
            for itk_child = BF_childidx(npkk,itk)
                len = length(skpre{xxparidx,itk_child});
                GTolcell{ell}{itx,itk_child} = UR(:, offset+(1:len));
                offset = offset+len;
            end
        end
    end
end

%---------------------------------------------------------------
%   U construction

npxx = npx*2^levels;
npkk = npk/2^levels;
Ucell = cell(prod(npxx),prod(npkk));
for itx = 1:prod(npxx)
    xxval = xxset{levels+1}{itx};
    for itk = 1:prod(npkk)
        Ucell{itx,itk} = fun(xxval,sk{itx,itk});
    end
end

%---------------------------------------------------------------
%   H construction

HTolcell = cell(levels,1);
for ell = 1:levels
    npkk = npk*2^ell;
    npxx = npx/2^ell;
    sxpre = sx;
    sx = cell(prod(npkk),prod(npxx));

    HTolcell{ell} = cell(prod(npkk),prod(npxx)*2^Dim);
    for itk = 1:prod(npkk)
        kkval = kkset{ell+1}{itk};
        kkparidx = BF_parentidx(npkk,itk);
        for itx = 1:prod(npxx)
            cursx = [];
            for itx_child = BF_childidx(npxx,itx)
                cursx = [cursx; sxpre{kkparidx,itx_child}];
            end
            [ CU, sx{itk,itx} ] = CURL(fun, cursx, kkval, r, tol);
            offset = 0;
            for itx_child = BF_childidx(npxx,itx)
                len = length(sxpre{kkparidx,itx_child});
                HTolcell{ell}{itk,itx_child} = CU(offset+(1:len),:);
                offset = offset+len;
            end
        end
    end
end

%---------------------------------------------------------------
%   V construction

npxx = npx/2^levels;
npkk = npk*2^levels;
Vcell = cell(prod(npkk),prod(npxx));
for itk = 1:prod(npkk)
    kkval = kkset{levels+1}{itk};
    for itx = 1:prod(npxx)
        Vcell{itk,itx} = fun(sx{itk,itx},kkval);
    end
end

%==============================================================
% Fast butterfly factorization compression
% The compression is split into two parts: outward compression
% and inward compression

nnz_init = BF_nnz(levels,Mcell,GTolcell,Ucell,HTolcell,Vcell);

[Mcell,GTolcell,Ucell,HTolcell,Vcell] = ...
    BF_inComp(levels,npx,npk,tol,...
    Mcell,GTolcell,Ucell,HTolcell,Vcell);
[Mcell,GTolcell,Ucell,HTolcell,Vcell] = ...
    BF_outComp(levels,npx,npk,tol,...
    Mcell,GTolcell,Ucell,HTolcell,Vcell);

nnz_comp = BF_nnz(levels,Mcell,GTolcell,Ucell,HTolcell,Vcell);

Rcomp = nnz_init/nnz_comp;

%==============================================================
% Sparse matrix assembling
%  K = U Gl ... G2 G1 M H1 H2 ... Hl V
% where U, Gi, M, Hi, and V are sparse matrices

%--------------------------------------------------------------
%   M assembling

totalel = 0;
for itx = 1:prod(npx)
    for itk = 1:prod(npk)
        totalel = totalel + numel(Mcell{itx,itk});
    end
end

totalH = 0;
offsetH = zeros(prod(npx),prod(npk));
for itx = 1:prod(npx)
    for itk = 1:prod(npk)
        offsetH(itx,itk) = totalH;
        totalH = totalH + size(Mcell{itx,itk},1);
    end
end

totalW = 0;
offsetW = zeros(prod(npk),prod(npx));
for itk = 1:prod(npk)
    for itx = 1:prod(npx)
        offsetW(itk,itx) = totalW;
        totalW = totalW + size(Mcell{itx,itk},2);
    end
end

XT = zeros(totalel,1);
YT = zeros(totalel,1);
ST = zeros(totalel,1);
offset = 0;
for itx = 1:prod(npx)
    for itk = 1:prod(npk)
        Mat = Mcell{itx,itk};
        [MatH,MatW] = size(Mat);
        X = (offsetH(itx,itk)+(1:MatH))'*ones(1,MatW);
        Y = ones(MatH,1)*(offsetW(itk,itx)+(1:MatW));
        idx = offset+(1:MatH*MatW);
        XT(idx) = X(:);
        YT(idx) = Y(:);
        ST(idx) = Mat(:);
        if(~isempty(idx))
            offset = idx(end);
        end
    end
end
Factor.M = sparse(XT,YT,ST,totalH,totalW);


%---------------------------------------------------------------
%   G assembling

Factor.GTol = cell(levels,1);
for ell = 1:levels
    npxx = npx*2^ell;
    npkk = npk/2^ell;
    npxx_par = npx*2^(ell-1);
    npkk_child = npk/2^(ell-1);
    totalel = 0;
    for itx = 1:prod(npxx)
        for itk_child = 1:prod(npkk_child)
            totalel = totalel + numel(GTolcell{ell}{itx,itk_child});
        end
    end
    
    totalH = 0;
    offsetH = zeros(prod(npxx),prod(npkk_child));
    itkchildcell = zeros(prod(npkk),2^Dim);
    for itk = 1:prod(npkk)
        for it_child = 1:2^Dim
            itkchildcell(itk,it_child) = ...
                BF_vec2idx(npkk_child,(BF_idx2vec(npkk,itk)-1)*2 ...
                +BF_idx2vec(2*ones(1,Dim),it_child));
        end
    end
    for itx = 1:prod(npxx)
        for itk = 1:prod(npkk)
            for it_child = 1:2^Dim
                itk_child = itkchildcell(itk,it_child);
                offsetH(itx,itk_child) = totalH;
            end
            totalH = totalH + size(GTolcell{ell}{itx,itk_child},1);
        end
    end
    
    totalW = 0;
    offsetW = zeros(prod(npxx),prod(npkk_child));
    itxcell = zeros(prod(npxx_par),2^Dim);
    for itx_par = 1:prod(npxx_par)
        for it_child = 1:2^Dim
            itxcell(itx_par,it_child) = ...
                BF_vec2idx(npxx,(BF_idx2vec(npxx_par,itx_par)-1)*2 ...
                +BF_idx2vec(2*ones(1,Dim),it_child));
        end
    end
    for itx_par = 1:prod(npxx_par)
        for itk_child = 1:prod(npkk_child)
            for it_child = 1:2^Dim
                itx = itxcell(itx_par,it_child);
                offsetW(itx,itk_child) = totalW;
            end
            totalW = totalW + size(GTolcell{ell}{itx,itk_child},2);
        end
    end
    
    XT = zeros(totalel,1);
    YT = zeros(totalel,1);
    ST = zeros(totalel,1);
    offset = 0;
    for itx = 1:prod(npxx)
        for itk_child = 1:prod(npkk_child)
            Mat = GTolcell{ell}{itx,itk_child};
            [MatH,MatW] = size(Mat);
            X = (offsetH(itx,itk_child)+(1:MatH))'*ones(1,MatW);
            Y = ones(MatH,1)*(offsetW(itx,itk_child)+(1:MatW));
            idx = offset+(1:MatH*MatW);
            XT(idx) = X(:);
            YT(idx) = Y(:);
            ST(idx) = Mat(:);
            if(~isempty(idx))
                offset = idx(end);
            end
        end
    end
    Factor.GTol{ell} = sparse(XT,YT,ST,totalH,totalW);
end

%---------------------------------------------------------------
%   U assembling

npxx = npx*2^levels;
npkk = npk/2^levels;
totalel = 0;
for itx = 1:prod(npxx)
    for itk = 1:prod(npkk)
        totalel = totalel + numel(Ucell{itx,itk});
    end
end

totalH = size(xx,1);

totalW = 0;
offsetW = zeros(prod(npxx),prod(npkk));
for itx = 1:prod(npxx)
    for itk = 1:prod(npkk)
        offsetW(itx,itk) = totalW;
        totalW = totalW + size(Ucell{itx,itk},2);
    end
end

XT = zeros(totalel,1);
YT = zeros(totalel,1);
ST = zeros(totalel,1);
offset = 0;
for itx = 1:prod(npxx)
    for itk = 1:prod(npkk)
        Mat = Ucell{itx,itk};
        [MatH,MatW] = size(Mat);
        X = (xxidx(xxIA(itx):xxIA(itx+1)-1))*ones(1,MatW);
        Y = ones(MatH,1)*(offsetW(itx,itk)+(1:MatW));
        idx = offset+(1:MatH*MatW);
        XT(idx) = X(:);
        YT(idx) = Y(:);
        ST(idx) = Mat(:);
        if(~isempty(idx))
            offset = idx(end);
        end
    end
end
Factor.U = sparse(XT,YT,ST,totalH,totalW);

%---------------------------------------------------------------
%   H assembling

Factor.HTol = cell(levels,1);
for ell = 1:levels
    npkk = npk*2^ell;
    npxx = npx/2^ell;
    npkk_par = npk*2^(ell-1);
    npxx_child = npx/2^(ell-1);
    totalel = 0;
    for itk = 1:prod(npkk)
        for itx_child = 1:prod(npxx_child)
            totalel = totalel + numel(HTolcell{ell}{itk,itx_child});
        end
    end
    
    totalW = 0;
    offsetW = zeros(prod(npkk),prod(npxx_child));
    itxchildcell = zeros(prod(npxx),2^Dim);
    for itx = 1:prod(npxx)
        for it_child = 1:2^Dim
            itxchildcell(itx,it_child) = ...
                BF_vec2idx(npxx_child,(BF_idx2vec(npxx,itx)-1)*2 ...
                +BF_idx2vec(2*ones(1,Dim),it_child));
        end
    end
    for itk = 1:prod(npkk)
        for itx = 1:prod(npxx)
            for it_child = 1:2^Dim
                itx_child = itxchildcell(itx,it_child);
                offsetW(itk,itx_child) = totalW;
            end
            totalW = totalW + size(HTolcell{ell}{itk,itx_child},2);
        end
    end
    
    totalH = 0;
    offsetH = zeros(prod(npkk),prod(npxx_child));
    itkcell = zeros(prod(npkk_par),2^Dim);
    for itk_par = 1:prod(npkk_par)
        for it_child = 1:2^Dim
            itkcell(itk_par,it_child) = ...
                BF_vec2idx(npkk,(BF_idx2vec(npkk_par,itk_par)-1)*2 ...
                +BF_idx2vec(2*ones(1,Dim),it_child));
        end
    end
    for itk_par = 1:prod(npkk_par)
        for itx_child = 1:prod(npxx_child)
            for it_child = 1:2^Dim
                itk = itkcell(itk_par,it_child);
                offsetH(itk,itx_child) = totalH;
            end
            totalH = totalH + size(HTolcell{ell}{itk,itx_child},1);
        end
    end
    
    XT = zeros(totalel,1);
    YT = zeros(totalel,1);
    ST = zeros(totalel,1);
    offset = 0;
    for itk = 1:prod(npkk)
        for itx_child = 1:prod(npxx_child)
            Mat = HTolcell{ell}{itk,itx_child};
            [MatH,MatW] = size(Mat);
            X = (offsetH(itk,itx_child)+(1:MatH))'*ones(1,MatW);
            Y = ones(MatH,1)*(offsetW(itk,itx_child)+(1:MatW));
            idx = offset+(1:MatH*MatW);
            XT(idx) = X(:);
            YT(idx) = Y(:);
            ST(idx) = Mat(:);
            if(~isempty(idx))
                offset = idx(end);
            end
        end
    end
    Factor.HTol{ell} = sparse(XT,YT,ST,totalH,totalW);
end

%---------------------------------------------------------------
%   V assembling

npkk = npk*2^levels;
npxx = npx/2^levels;
totalel = 0;
for itk = 1:prod(npkk)
    for itx = 1:prod(npxx)
        totalel = totalel + numel(Vcell{itk,itx});
    end
end

totalH = 0;
offsetH = zeros(prod(npkk),prod(npxx));
for itk = 1:prod(npkk)
    for itx = 1:prod(npxx)
        offsetH(itk,itx) = totalH;
        totalH = totalH + size(Vcell{itk,itx},1);
    end
end

totalW = size(kk,1);

XT = zeros(totalel,1);
YT = zeros(totalel,1);
ST = zeros(totalel,1);
offset = 0;
for itk = 1:prod(npkk)
    for itx = 1:prod(npxx)
        Mat = Vcell{itk,itx};
        [MatH,MatW] = size(Mat);
        X = (offsetH(itk,itx)+(1:MatH))'*ones(1,MatW);
        Y = ones(MatH,1)*(kkidx(kkIA(itk):kkIA(itk+1)-1))';
        idx = offset+(1:MatH*MatW);
        XT(idx) = X(:);
        YT(idx) = Y(:);
        ST(idx) = Mat(:);
        if(~isempty(idx))
            offset = idx(end);
        end
    end
end
Factor.V = sparse(XT,YT,ST,totalH,totalW);

end
