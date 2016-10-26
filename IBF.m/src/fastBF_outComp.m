function [Mcell,GTolcell,Ucell,HTolcell,Vcell] = ...
    fastBF_outComp(levels,npx,npk,tol,Mcell,GTolcell,Ucell,HTolcell,Vcell)

Dim = size(npx,2);

%---------------------------------------------------------------
%   M compression

Lcell = cell(prod(npx),prod(npk));
Rcell = cell(prod(npx),prod(npk));

for itx = 1:prod(npx)
    for itk = 1:prod(npk)
        [U,S,V] = svdtrunc(Mcell{itx,itk},tol);
        Lcell{itx,itk} = U*sqrt(S);
        Rcell{itk,itx} = sqrt(S)*V';
        Mcell{itx,itk} = eye(size(S));
    end
end

%---------------------------------------------------------------
%   G compression

for ell = 1:levels
    npxx = npx*2^ell;
    npkk = npk/2^ell;
    npxx_par = npx*2^(ell-1);
    npkk_child = npk/2^(ell-1);
    for itx = 1:prod(npxx)
        itx_par = vec2idx(npxx_par,floor((idx2vec(npxx,itx)+1)/2));
        for itk_child = 1:prod(npkk_child)
            GTolcell{ell}{itx,itk_child} = ...
                GTolcell{ell}{itx,itk_child}*Lcell{itx_par,itk_child};
        end
    end
    
    Lcell = cell(prod(npxx),prod(npkk));
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
            tmpMat = [];
            for it_child = 1:2^Dim
                itk_child = itkchildcell(itk,it_child);
                tmpMat = [ tmpMat GTolcell{ell}{itx,itk_child} ];
            end
            [U,S,V] = svdtrunc(tmpMat,tol);
            Lcell{itx,itk} = U*S;
            offset = 0;
            for it_child = 1:2^Dim
                itk_child = itkchildcell(itk,it_child);
                wid = size(GTolcell{ell}{itx,itk_child},2);
                GTolcell{ell}{itx,itk_child} = V(offset+(1:wid),:)';
                offset = offset+wid;
            end
        end
    end
end


%---------------------------------------------------------------
%   U compression

npxx = npx*2^levels;
npkk = npk/2^levels;
for itx = 1:prod(npxx)
    for itk = 1:prod(npkk)
        Ucell{itx,itk} = Ucell{itx,itk}*Lcell{itx,itk};
    end
end


%---------------------------------------------------------------
%   H compression

for ell = 1:levels
    npkk = npk*2^ell;
    npxx = npx/2^ell;
    npkk_par = npk*2^(ell-1);
    npxx_child = npx/2^(ell-1);
    for itk = 1:prod(npkk)
        itk_par = vec2idx(npkk_par,floor((idx2vec(npkk,itk)+1)/2));
        for itx_child = 1:prod(npxx_child)
            HTolcell{ell}{itk,itx_child} = ...
                Rcell{itk_par,itx_child}*HTolcell{ell}{itk,itx_child};
        end
    end
    
    Rcell = cell(prod(npkk),prod(npxx));
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
            tmpMat = [];
            for it_child = 1:2^Dim
                itx_child = itxchildcell(itx,it_child);
                tmpMat = [ tmpMat; HTolcell{ell}{itk,itx_child} ];
            end
            [U,S,V] = svdtrunc(tmpMat,tol);
            Rcell{itk,itx} = S*V';
            offset = 0;
            for it_child = 1:2^Dim
                itx_child = itxchildcell(itx,it_child);
                hgh = size(HTolcell{ell}{itk,itx_child},1);
                HTolcell{ell}{itk,itx_child} = U(offset+(1:hgh),:);
                offset = offset+hgh;
            end
        end
    end
end


%---------------------------------------------------------------
%   V compression

npkk = npk*2^levels;
npxx = npx/2^levels;
for itk = 1:prod(npkk)
    for itx = 1:prod(npxx)
        Vcell{itk,itx} = Rcell{itk,itx}*Vcell{itk,itx};
    end
end

end