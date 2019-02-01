function B = HSSBF_sp_embedding(BF,diagFac)
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
%
% Reference: A Hierarchical Butterfly LU Preconditioner for Two-Dimensional
% Electromagnetic Scattering Problems Involving Open Surfaces, preprint,
% 2019
% Copyright reserved.
% By Haizhao Yang, 2018, National University of Singapore

if nargin < 2
    isD = 0;
else
    isD = length(diagFac);
    lenD = isD;
end

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

end
