function [Mcell,GTolcell,Ucell,HTolcell,Vcell] = ...
    fastBF_inComp(levels,npx,npk,tol,Mcell,GTolcell,Ucell,HTolcell,Vcell)

Dim = size(npx,2);


%---------------------------------------------------------------
%   U compression

npxx = npx*2^levels;
npkk = npk/2^levels;
Lcell = cell(prod(npxx),prod(npkk));
for itx = 1:prod(npxx)
    for itk = 1:prod(npkk)
        [U,S,V] = svdtrunc(Ucell{itx,itk},tol);
        Ucell{itx,itk} = U;
        Lcell{itx,itk} = S*V';
    end
end


%---------------------------------------------------------------
%   G compression

for ell = levels:-1:1
    npxx = npx*2^ell;
    npkk = npk/2^ell;
    npxx_par = npx*2^(ell-1);
    npkk_child = npk/2^(ell-1);
    for itx = 1:prod(npxx)
        for itk = 1:prod(npkk)
            for it_child = 1:2^Dim
                itk_child = vec2idx(npkk_child,(idx2vec(npkk,itk)-1)*2 ...
                    +idx2vec(2*ones(1,Dim),it_child));
                GTolcell{ell}{itx,itk_child} = ...
                    Lcell{itx,itk}*GTolcell{ell}{itx,itk_child};
            end
        end
    end
    
    Lcell = cell(prod(npxx_par),prod(npkk_child));
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
            tmpMat = [];
            for it_child = 1:2^Dim
                itx = itxcell(itx_par,it_child);
                tmpMat = [tmpMat; GTolcell{ell}{itx,itk_child}];
            end
            [U,S,V] = svdtrunc(tmpMat,tol);
            Lcell{itx_par,itk_child} = S*V';
            offset = 0;
            for it_child = 1:2^Dim
                itx = itxcell(itx_par,it_child);
                hgh = size(GTolcell{ell}{itx,itk_child},1);
                GTolcell{ell}{itx,itk_child} = U(offset+(1:hgh),:);
                offset = offset+hgh;
            end
        end
    end
end


%---------------------------------------------------------------
%   V compression

npkk = npk*2^levels;
npxx = npx/2^levels;
Rcell = cell(prod(npkk),prod(npxx));
for itk = 1:prod(npkk)
    for itx = 1:prod(npxx)
        [U,S,V] = svdtrunc(Vcell{itk,itx},tol);
        Vcell{itk,itx} = V';
        Rcell{itk,itx} = U*S;
    end
end


%---------------------------------------------------------------
%   H compression

for ell = levels:-1:1
    npkk = npk*2^ell;
    npxx = npx/2^ell;
    npkk_par = npk*2^(ell-1);
    npxx_child = npx/2^(ell-1);
    for itk = 1:prod(npkk)
        for itx = 1:prod(npxx)
            for it_child = 1:2^Dim
                itx_child = vec2idx(npxx_child,(idx2vec(npxx,itx)-1)*2 ...
                    +idx2vec(2*ones(1,Dim),it_child));
                HTolcell{ell}{itk,itx_child} = ...
                    HTolcell{ell}{itk,itx_child}*Rcell{itk,itx};
            end
        end
    end
    
    Rcell = cell(prod(npkk_par),prod(npxx_child));
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
            tmpMat = [];
            for it_child = 1:2^Dim
                itk = itkcell(itk_par,it_child);
                tmpMat = [tmpMat HTolcell{ell}{itk,itx_child}];
            end
            [U,S,V] = svdtrunc(tmpMat,tol);
            Rcell{itk_par,itx_child} = U*S;
            offset = 0;
            for it_child = 1:2^Dim
                itk = itkcell(itk_par,it_child);
                wid = size(HTolcell{ell}{itk,itx_child},2);
                HTolcell{ell}{itk,itx_child} = V(offset+(1:wid),:)';
                offset = offset+wid;
            end
        end
    end
end


%---------------------------------------------------------------
%   M compression

for itx = 1:prod(npx)
    for itk = 1:prod(npk)
        Mcell{itx,itk} = Lcell{itx,itk}*Rcell{itk,itx};
    end
end

end