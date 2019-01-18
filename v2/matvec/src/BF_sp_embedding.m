function B = BF_sp_embedding(BF,diagFac)
% This code embeds the butterfly factorization (stored in BF) into a
% large sparse matrix. In the presence of diagonal factors, put it in the
% diagonal part of B.
%
% Reference: Equation (2.4) in
% A FAST SEMI-DIRECT LEAST SQUARES ALGORITHM FOR HIERARCHICALLY BLOCK
% SEPARABLE MATRICES.
% KENNETH L. HO AND LESLIE GREENGARD, SIMAX, 2014.
% BF stores matrix factors L^{i}, R^{i}, and D^{0}
% diagFac stroes matrix factors D^{i} for i \neq 0
%
% Copyright 2018 Haizhao Yang
if nargin < 2
    isD = 0;
else
    isD = length(diagFac);
    lenD = isD;
end

if iscell(BF)
    fn = length(BF);
    hn = (length(BF)-1)/2;
    gn = hn;
    
    len = 0;
    for cnt = 1:fn
        len = len + size(BF{cnt},2);
    end
    B = sparse(len,len);
    
    cedU = 0;
    redV = 0;
    for cnt = 1:gn
        cedUo = cedU;
        redVo = redV;
        U = BF{cnt};
        V = BF{fn-cnt+1};
        rstU = 1+redVo;
        redU = size(U,1)+redVo;
        cstU = 1+size(V,2)+cedUo;
        cedU = size(U,2)+size(V,2)+cedUo;
        B(rstU:redU,cstU:cedU) = U;
        rstV = 1+size(U,1)+redVo;
        redV = size(V,1)+size(U,1)+redVo;
        cstV = 1+cedUo;
        cedV = size(V,2)+cedUo;
        B(rstV:redV,cstV:cedV) = V;
        sz = size(BF{fn-cnt+1-1},2);
        B(rstV:redV,(cedU+1):(cedU+sz)) = -speye(size(V,1),sz);
        sz = size(BF{cnt+1},1);
        B((redV+1):(redV+sz),cstU:cedU) = -speye(sz,size(U,2));
        
        if isD > 0
            isD = isD - 1;
            B(rstU:redU,cstV:cedV) = diagFac{lenD-isD};
        end
    end
    
    if gn == 0
        redV = 0;
        cedU = 0;
    end
    S = BF{gn+1};
    B((redV+1):(redV+size(S,1)),(cedU+1):(cedU+size(S,2))) = S;
else
    
    len = 0;
    len = len + size(BF.V,2);
    for i=length(BF.HTol):-1:1
        len = len + size(BF.HTol{i},2);
    end
    len = len + size(BF.M,2);
    for i=1:length(BF.GTol)
        len = len + size(BF.GTol{i},2);
    end
    len = len + size(BF.U,2);
    B = sparse(len,len);
    
    
    hn = length(BF.HTol);
    gn = length(BF.GTol);
    
    rstU = 1;
    redU = size(BF.U,1);
    cstU = 1+size(BF.V,2);
    cedU = size(BF.U,2)+size(BF.V,2);
    B(rstU:redU,cstU:cedU) = BF.U;
    rstV = 1+size(BF.U,1);
    redV = size(BF.V,1)+size(BF.U,1);
    cstV = 1;
    cedV = size(BF.V,2);
    B(rstV:redV,cstV:cedV) = BF.V;
    if isD > 0
        isD = isD - 1;
        B(rstU:redU,cstV:cedV) = diagFac{lenD-isD};
    end
    if hn > 0
        sz = size(BF.HTol{hn},2);
        B(rstV:redV,(cedU+1):(cedU+sz)) = -speye(size(BF.V,1),sz);
    else
        sz = size(BF.M,2);
        B(rstV:redV,(cedU+1):(cedU+sz)) = -speye(size(BF.V,1),sz);
    end
    if gn > 0
        sz = size(BF.GTol{gn},1);
        B((redV+1):(redV+sz),cstU:cedU) = -speye(sz,size(BF.U,2));
    else
        sz = size(BF.M,1);
        B((redV+1):(redV+sz),cstU:cedU) = -speye(sz,size(BF.U,2));
    end
    
    if gn > 0
        cedUo = cedU;
        redVo = redV;
        for cnt = gn:-1:2
            U = BF.GTol{cnt};
            V = BF.HTol{cnt};
            rstU = 1+redVo;
            redU = size(U,1)+redVo;
            cstU = 1+size(V,2)+cedUo;
            cedU = size(U,2)+size(V,2)+cedUo;
            B(rstU:redU,cstU:cedU) = U;
            rstV = 1+size(U,1)+redVo;
            redV = size(V,1)+size(U,1)+redVo;
            cstV = 1+cedUo;
            cedV = size(V,2)+cedUo;
            B(rstV:redV,cstV:cedV) = V;
            sz = size(BF.HTol{cnt-1},2);
            B(rstV:redV,(cedU+1):(cedU+sz)) = -speye(size(V,1),sz);
            sz = size(BF.GTol{cnt-1},1);
            B((redV+1):(redV+sz),cstU:cedU) = -speye(sz,size(U,2));
            
            if isD > 0
                isD = isD - 1;
                B(rstU:redU,cstV:cedV) = diagFac{lenD-isD};
            end
        end
        
        cedUo = cedU;
        redVo = redV;
        cnt = 1;
        U = BF.GTol{cnt};
        V = BF.HTol{cnt};
        rstU = 1+redVo;
        redU = size(U,1)+redVo;
        cstU = 1+size(V,2)+cedUo;
        cedU = size(U,2)+size(V,2)+cedUo;
        B(rstU:redU,cstU:cedU) = U;
        rstV = 1+size(U,1)+redVo;
        redV = size(V,1)+size(U,1)+redVo;
        cstV = 1+cedUo;
        cedV = size(V,2)+cedUo;
        B(rstV:redV,cstV:cedV) = V;
        sz = size(BF.M,2);
        B(rstV:redV,(cedU+1):(cedU+sz)) = -speye(size(V,1),sz);
        sz = size(BF.M,1);
        B((redV+1):(redV+sz),cstU:cedU) = -speye(sz,size(U,2));
        
        if isD > 0
            isD = isD - 1;
            B(rstU:redU,cstV:cedV) = diagFac{lenD-isD};
        end
    end
    S = BF.M;
    B((redV+1):(redV+size(S,1)),(cedU+1):(cedU+size(S,2))) = S;
end
end
