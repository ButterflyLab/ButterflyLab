function [Factor,Rcomp] = fastBF(fun_org,xx,kk,NG,tol,tag)

if nargin == 5
    tag = 'Regular';
end

[fun,xx,kk,xbox,kbox,npx,npk] = fbf_prepbox(fun_org,xx,kk,tag);

grid = Chey_grid(NG);
Dim = size(xx,2);
NGD = NG^Dim;

levels = max(0,min(floor(log2([npx npk])-3)));

xxboxidx = zeros(size(xx));
npxx = npx*2^levels;
for i = 1:Dim
    edges = linspace(xbox(1,i),xbox(2,i),npxx(i)+1);
    [~,xxboxidx(:,i)] = histc(xx(:,i),edges);
end
[xxboxidx,xxidx] = sortrows(xxboxidx,Dim:-1:1);
[xxC,xxIA,~] = unique(xxboxidx,'rows','stable');
xxC = [xxC;zeros(1,size(xxC,2))];
xxIA = [xxIA;size(xx,1)+1];
for itx = 1:prod(npxx)
    x = idx2vec(npxx,itx);
    if any(x ~= xxC(itx,:))
        xxC = [xxC(1:itx-1,:);x;xxC(itx:end,:)];
        xxIA = [xxIA(1:itx-1);xxIA(itx);xxIA(itx:end)];
    end
end

kkboxidx = zeros(size(kk));
npkk = npk*2^levels;
for i = 1:Dim
    edges = linspace(kbox(1,i),kbox(2,i),npkk(i)+1);
    [~,kkboxidx(:,i)] = histc(kk(:,i),edges);
end
[kkboxidx,kkidx] = sortrows(kkboxidx,Dim:-1:1);
[kkC,kkIA,~] = unique(kkboxidx,'rows','stable');
kkC = [kkC;zeros(1,size(kkC,2))];
kkIA = [kkIA;size(kk,1)+1];
for itk = 1:prod(npkk)
    k = idx2vec(npkk,itk);
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
    x = idx2vec(npx,itx);
    xgridcell{itx} = fbf_grid(x,npx,grid,xbox);
end
for itk = 1:prod(npk)
    k = idx2vec(npk,itk);
    kgridcell{itk} = fbf_grid(k,npk,grid,kbox);
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
        [xgrid,~] = fbf_grid(idx2vec(2*ones(1,Dim), ...
            it_child),npxx,grid,xbox);
        [~,xparLgrid] = fbf_grid(ones(1,Dim),npxx/2,grid,xbox);
        LagrangeMatCell{it_child} = fbf_Lagrange(xparLgrid,xgrid).';
    end
    
    GTolcell{ell} = cell(prod(npxx),prod(npkk_child));
    kchildcell = cell(prod(npkk_child),1);
    for itk_child = 1:prod(npkk_child)
        kchildcell{itk_child} = idx2vec(npkk_child,itk_child);
    end
    for itx = 1:prod(npxx)
        x = idx2vec(npxx,itx);
        [xgrid,~] = fbf_grid(x,npxx,grid,xbox);
        [xpargrid,~] = fbf_grid(floor((x-1)/2)+1,npxx/2,grid,xbox);
        it_child = vec2idx(2*ones(1,Dim),x-2*floor((x-1)/2));
        LagrangeMat = LagrangeMatCell{it_child};
        for itk_child = 1:prod(npkk_child)
            k_child = kchildcell{itk_child};
            kcen_child = kbox(1,:) ...
                +(k_child-1/2).*(kbox(2,:)-kbox(1,:))./npkk_child;
            GTolcell{ell}{itx,itk_child} = ...
                sparse(1:NGD,1:NGD,fun(xgrid,kcen_child))*(LagrangeMat*...
                sparse(1:NGD,1:NGD,1./fun(xpargrid,kcen_child)));
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
    kcell{itk} = idx2vec(npkk,itk);
end
for itx = 1:prod(npxx)
    x = idx2vec(npxx,itx);
    [xgrid,xLgrid] = fbf_grid(x,npxx,grid,xbox);
    xxsub = xx(xxidx(xxIA(itx):xxIA(itx+1)-1),:);
    LagrangeMat = fbf_Lagrange(xLgrid,xxsub).';
    for itk = 1:prod(npkk)
        k = kcell{itk};
        kcen = kbox(1,:)+(k-1/2).*(kbox(2,:)-kbox(1,:))./npkk;
        Ucell{itx,itk} = ...
            sparse(1:size(xxsub,1),1:size(xxsub,1),fun(xxsub,kcen))...
            *(LagrangeMat*(sparse(1:NGD,1:NGD,1./fun(xgrid,kcen))));
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
        [kgrid,~] = fbf_grid(idx2vec(2*ones(1,Dim),it_child),...
            npkk,grid,kbox);
        [~,kparLgrid] = fbf_grid(ones(1,Dim),npkk/2,grid,kbox);
        LagrangeMatCell{it_child} = fbf_Lagrange(kparLgrid,kgrid);
    end
    
    HTolcell{ell} = cell(prod(npkk),prod(npxx_child));
    xchildcell = cell(prod(npxx_child),1);
    for itx_child = 1:prod(npxx_child)
        xchildcell{itx_child} = idx2vec(npxx_child,itx_child);
    end
    for itk = 1:prod(npkk)
        k = idx2vec(npkk,itk);
        kgrid = fbf_grid(k,npkk,grid,kbox);
        [kpargrid,~] = fbf_grid(floor((k-1)/2)+1,npkk/2,grid,kbox);
        it_child = vec2idx(2*ones(1,Dim),k-2*floor((k-1)/2));
        LagrangeMat = LagrangeMatCell{it_child};
        for itx_child = 1:prod(npxx_child)
            x_child = xchildcell{itx_child};
            xcen_child = xbox(1,:)+(x_child-1/2).*(xbox(2,:) ...
                -xbox(1,:))./npxx_child;
            HTolcell{ell}{itk,itx_child} = ...
                sparse(1:NGD,1:NGD,1./fun(xcen_child,kpargrid)) ...
                * LagrangeMat * sparse(1:NGD,1:NGD,fun(xcen_child,kgrid));
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
    xcell{itx} = idx2vec(npxx,itx);
end
for itk = 1:prod(npkk)
    k = idx2vec(npkk,itk);
    [kgrid,kLgrid] = fbf_grid(k,npkk,grid,kbox);
    kksub = kk(kkidx(kkIA(itk):kkIA(itk+1)-1),:);
    LagrangeMat = fbf_Lagrange(kLgrid,kksub);
    for itx = 1:prod(npxx)
        x = xcell{itx};
        xcen = xbox(1,:)+(x-1/2).*(xbox(2,:)-xbox(1,:))./npxx;
        Vcell{itk,itx} = ...
            sparse(1:NGD,1:NGD,1./fun(xcen,kgrid))*LagrangeMat*...
            sparse(1:size(kksub,1),1:size(kksub,1),fun(xcen,kksub));
    end
end

%==============================================================
% Fast butterfly factorization compression
% The compression is split into two parts: outward compression
% and inward compression

nnz_init = fbf_nnz(levels,Mcell,GTolcell,Ucell,HTolcell,Vcell);

[Mcell,GTolcell,Ucell,HTolcell,Vcell] = ...
    fastBF_outComp(levels,npx,npk,tol,...
    Mcell,GTolcell,Ucell,HTolcell,Vcell);
[Mcell,GTolcell,Ucell,HTolcell,Vcell] = ...
    fastBF_inComp(levels,npx,npk,tol,...
    Mcell,GTolcell,Ucell,HTolcell,Vcell);

nnz_comp = fbf_nnz(levels,Mcell,GTolcell,Ucell,HTolcell,Vcell);

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
                vec2idx(npkk_child,(idx2vec(npkk,itk)-1)*2 ...
                +idx2vec(2*ones(1,Dim),it_child));
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
                vec2idx(npxx,(idx2vec(npxx_par,itx_par)-1)*2 ...
                +idx2vec(2*ones(1,Dim),it_child));
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
                vec2idx(npxx_child,(idx2vec(npxx,itx)-1)*2 ...
                +idx2vec(2*ones(1,Dim),it_child));
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
                vec2idx(npkk,(idx2vec(npkk_par,itk_par)-1)*2 ...
                +idx2vec(2*ones(1,Dim),it_child));
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
