function [Factor,Rcomp] = IBF_uniform(fun,xx,kk,NG,tol,disSet,isRecomp,tag)
% O(N log N) operation and memory complexity.
% Reference: H. Yang, A Unified Framework for Oscillatory Integral 
% Transform: When to use NUFFT or Butterfly factorization? preprint, 2018.
%
% Inputs:
% fun is a function handle such that fun(xx,kk) is the complementary
% low-rank matrix to be compressed.
% xx is a column vector and must be (1:numel(xx))'.
% kk is a column vector and must be (1:numel(kk))'.
% NG is rank of low-rank matrices
% disSet - the set of discontinuous points of the phase function of fun.
% isRecomp - 1: recompress the butterfly factorization using the sweeping
% compression technique; 0: no recompression.
%
% By Haizhao Yang, 2018

if nargin < 6, disSet = []; end
if nargin < 7, isRecomp = 1; end
if nargin < 8
    tag = 'Regular';
end


if size(xx,1)<size(xx,2), error('xx should be a tall skinny matrix'); end
if size(kk,1)<size(kk,2), error('kk should be a tall skinny matrix'); end

Dim = 1;
grid = BF_Chey_grid(NG);
[xbox,kbox,npx,npk] = BF_prepbox_int(numel(xx),numel(kk),tag);

% make sure that the smallest dimension of the low-rank submatrices is
% larger than max(8,NG)
levels = max(0,max(floor(log2([npx npk])-max(3,ceil(log2(NG)))))); 

npxx = npx*2^levels;
edges = linspace(xbox(1),xbox(2),npxx+1); % partition the range [xbox(1),xbox(2)] into npxx boxes
edges(end) = edges(end)+1;
[~,xxboxidx] = histc(xx(:),edges);
[xxboxidx,xxidx] = sortrows(xxboxidx,1);
[xxC,xxIA,~] = unique(xxboxidx,'rows','stable');
xxC = [xxC;zeros(1,size(xxC,2))];
xxIA = [xxIA;size(xx,1)+1]; % the start point of the i-th box
for itx = 1:npxx % loop for each box
    x = BF_idx2vec(npxx,itx);
    if any(x ~= xxC(itx,:))
        xxC = [xxC(1:itx-1,:);x;xxC(itx:end,:)];
        xxIA = [xxIA(1:itx-1);xxIA(itx);xxIA(itx:end)];
    end
end

npkk = npk*2^levels;
edges = linspace(kbox(1),kbox(2),npkk+1);
edges(end) = edges(end)+1;
[~,kkboxidx] = histc(kk(:),edges);
[kkboxidx,kkidx] = sortrows(kkboxidx,1);
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
xgridcell = cell(prod(npx),1);
kgridcell = cell(prod(npk),1);
for itx = 1:prod(npx)
    x = BF_idx2vec(npx,itx);
    xgridcell{itx} = BF_grid_int(grid,x,NG,npx,xbox,disSet);
end
for itk = 1:prod(npk)
    k = BF_idx2vec(npk,itk);
    kgridcell{itk} = BF_grid_int(grid,k,NG,npk,kbox,disSet);
end

for itx = 1:prod(npx)
    for itk = 1:prod(npk)
        xgrid = xgridcell{itx};
        kgrid = kgridcell{itk};
        Mcell{itx,itk} = fun(xgrid,kgrid);
    end
end

%---------------------------------------------------------------
%   G construction


GTolcell = cell(levels,1);
for ell = 1:levels
    npxx = npx*2^ell;
    npkk_child = npk/2^(ell-1);
    
    LagrangeMatCell = cell(2^Dim);
    for it_child = 1:2^Dim
        xgrid = BF_grid_int(grid,BF_idx2vec(2*ones(1,Dim), ...
            it_child),NG,npxx,xbox,disSet);
        xparLgrid = BF_grid_int(grid,ones(1,Dim),NG,npxx/2,xbox,disSet);
        LagrangeMatCell{it_child} = BF_Lagrange(xparLgrid(:),xgrid(:)).';
    end
    
    GTolcell{ell} = cell(prod(npxx),prod(npkk_child));
    kchildcell = cell(prod(npkk_child),1);
    for itk_child = 1:prod(npkk_child)
        kchildcell{itk_child} = BF_idx2vec(npkk_child,itk_child);
    end
    for itx = 1:prod(npxx)
        x = BF_idx2vec(npxx,itx);
        xgrid = BF_grid_int(grid,x,NG,npxx,xbox,disSet);
        xpargrid = BF_grid_int(grid,floor((x-1)/2)+1,NG,npxx/2,xbox,disSet);
        it_child = BF_vec2idx(2*ones(1,Dim),x-2*floor((x-1)/2));
        LagrangeMat = LagrangeMatCell{it_child};
        for itk_child = 1:prod(npkk_child)
            k_child = kchildcell{itk_child};
            kcen_child = round(kbox(1,:) ...
                +(k_child-1/2).*(kbox(2,:)-kbox(1,:))./npkk_child);
            GTolcell{ell}{itx,itk_child} = ...
                sparse(1:numel(xgrid),1:numel(xgrid),fun(xgrid,kcen_child))*(LagrangeMat*...
                sparse(1:numel(xpargrid),1:numel(xpargrid),1./fun(xpargrid,kcen_child)));
        end
    end
end

%---------------------------------------------------------------
%   U construction

npxx = npx*2^levels;
npkk = npk/2^levels;
Ucell = cell(prod(npxx),prod(npkk));
kcell = cell(prod(npkk),1);
for itk = 1:prod(npkk)
    kcell{itk} = BF_idx2vec(npkk,itk);
end
for itx = 1:prod(npxx)
    x = BF_idx2vec(npxx,itx);
    xgrid = BF_grid_int(grid,x,NG,npxx,xbox,disSet); xLgrid = xgrid;
    xxsub = xx(xxidx(xxIA(itx):xxIA(itx+1)-1));
    LagrangeMat = BF_Lagrange(xLgrid(:),xxsub(:)).';
    for itk = 1:prod(npkk)
        k = kcell{itk};
        kcen = round(kbox(1,:)+(k-1/2).*(kbox(2,:)-kbox(1,:))./npkk);
        Ucell{itx,itk} = ...
            sparse(1:numel(xxsub),1:numel(xxsub),fun(xxsub,kcen))...
            *(LagrangeMat*(sparse(1:numel(xgrid),1:numel(xgrid),1./fun(xgrid,kcen))));
    end
end

%---------------------------------------------------------------
%   H construction

HTolcell = cell(levels,1);
for ell = 1:levels
    npkk = npk*2^ell;
    npxx_child = npx/2^(ell-1);
    
    LagrangeMatCell = cell(2^Dim);
    for it_child = 1:2^Dim
        kgrid = BF_grid_int(grid,BF_idx2vec(2*ones(1,Dim),it_child),...
            NG,npkk,kbox,disSet);
        kparLgrid = BF_grid_int(grid,ones(1,Dim),NG,npkk/2,kbox,disSet);
        LagrangeMatCell{it_child} = BF_Lagrange(kparLgrid(:),kgrid(:));
    end
    
    HTolcell{ell} = cell(prod(npkk),prod(npxx_child));
    xchildcell = cell(prod(npxx_child),1);
    for itx_child = 1:prod(npxx_child)
        xchildcell{itx_child} = BF_idx2vec(npxx_child,itx_child);
    end
    for itk = 1:prod(npkk)
        k = BF_idx2vec(npkk,itk);
        kgrid = BF_grid_int(grid,k,NG,npkk,kbox,disSet);
        kpargrid = BF_grid_int(grid,floor((k-1)/2)+1,NG,npkk/2,kbox,disSet);
        it_child = BF_vec2idx(2*ones(1,Dim),k-2*floor((k-1)/2));
        LagrangeMat = LagrangeMatCell{it_child};
        for itx_child = 1:prod(npxx_child)
            x_child = xchildcell{itx_child};
            xcen_child = round(xbox(1,:)+(x_child-1/2).*(xbox(2,:) ...
                -xbox(1,:))./npxx_child);
            HTolcell{ell}{itk,itx_child} = ...
                sparse(1:numel(kpargrid),1:numel(kpargrid),1./fun(xcen_child,kpargrid)) ...
                * LagrangeMat * sparse(1:numel(kgrid),1:numel(kgrid),fun(xcen_child,kgrid));
        end
    end
end

%---------------------------------------------------------------
%   V construction

npkk = npk*2^levels;
npxx = npx/2^levels;
Vcell = cell(prod(npkk),prod(npxx));
xcell = cell(prod(npxx),1);
for itx = 1:prod(npxx)
    xcell{itx} = BF_idx2vec(npxx,itx);
end
for itk = 1:prod(npkk)
    k = BF_idx2vec(npkk,itk);
    kgrid = BF_grid_int(grid,k,NG,npkk,kbox,disSet); kLgrid = kgrid;
    kksub = kk(kkidx(kkIA(itk):kkIA(itk+1)-1));
    LagrangeMat = BF_Lagrange(kLgrid(:),kksub(:));
    for itx = 1:prod(npxx)
        x = xcell{itx};
        xcen = round(xbox(1,:)+(x-1/2).*(xbox(2,:)-xbox(1,:))./npxx);
        Vcell{itk,itx} = ...
            sparse(1:numel(kgrid),1:numel(kgrid),1./fun(xcen,kgrid))*LagrangeMat*...
            sparse(1:numel(kksub),1:numel(kksub),fun(xcen,kksub));
    end
end

%==============================================================
% Fast butterfly factorization compression
% The compression is split into two parts: outward compression
% and inward compression
if isRecomp
    nnz_init = BF_nnz(levels,Mcell,GTolcell,Ucell,HTolcell,Vcell);
    
    [Mcell,GTolcell,Ucell,HTolcell,Vcell] = ...
        BF_outComp(levels,npx,npk,tol,...
        Mcell,GTolcell,Ucell,HTolcell,Vcell);
    [Mcell,GTolcell,Ucell,HTolcell,Vcell] = ...
        BF_inComp(levels,npx,npk,tol,...
        Mcell,GTolcell,Ucell,HTolcell,Vcell);
    
    nnz_comp = BF_nnz(levels,Mcell,GTolcell,Ucell,HTolcell,Vcell);
    
    Rcomp = nnz_init/nnz_comp;
else
    Rcomp = 1;
end

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

sparse(XT,YT,ST,totalH,totalW);
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
