function [U,S,V,Ua,Sa,Va] = MVBF_Lowrank_MD(Nx,Np,fun,funt,tol,tR,mR,dim,thre,type,isCmpl,opt)
%
% Input:
% fun is the compressed butterfly matrix in a form of a function handle
% funt is the (conjugate) transpose of fun in a form of a function handle
% [Nx,Np] is the size of the matrix fun
% tol is the acccuracy tollerence
% tR is the test rank
% mR is the maximum rank, tR>=mR
% dim is the dimension of the problem
% thre is the threshold for discontinuity detection
% type, when funt is transpose, type = 0; otherwise, type = 1.
% isCmpl, 1, if the input function handle is complex; 0, otherwise
% opt is a structure storing the parameters for mode decomposition
%
% Output:
% Low-rank factorization s.t. \sum_k (Ua{k}*Sa{k}*Va{k}').*exp(1i*U{k}*S{k}*V{k}') \approx fun
%
% O(N log N) operation and memory complexity.
% Reference: H. Yang, A Unified Framework for Oscillatory Integral 
% Transform: When to use NUFFT or Butterfly factorization? preprint, 2018.

if nargin < 12 % set defaut parameters for mode decom[osition
    opt.numCom = 1; % TODO
    opt.VMDonly = 0;
end
if nargin < 11
    isCmpl = 0;
end
opt.isCmpl = isCmpl;
if nargin < 10
    type = 0;
end
if(Nx==0 || Np==0)
    U = zeros(Nx,0);
    S = zeros(0,0);
    V = zeros(Np,0);
    return;
end

% test discontinuity along row
% test discontinuity along column
% TODO: mode decomposition with discontinuity

disPos1 = [1];
disPos2 = [1];
Ridx = cell(opt.numCom,1); Cidx = cell(opt.numCom,1);
Ridxa = cell(opt.numCom,1); Cidxa = cell(opt.numCom,1);

if(tR<Np && tR<Nx)
    %rows
    rs = BF_RandSample(Nx,tR);
    rs = [disPos1,disPos1+1,disPos1+2,rs]; rs = unique(rs); rs = sort(rs);
    disPosSub1 = BF_find_entry(rs,disPos1);
    tmp = zeros(Np,numel(rs));
    tmp(rs+((1:numel(rs))-1)*Np) = 1;
    data = funt(tmp);
    data = BF_ModeDecom(data,opt);
    M2alla = abs(data);
    if type == 0
        data = data./abs(data);
    else
        data = conj(data)./abs(data);
    end
    M2all = imag(log(data));
    
    %columns
    cs = BF_RandSample(Np,tR);
    cs = [disPos2,disPos2+1,disPos2+2,cs]; cs = unique(cs); cs = sort(cs);
    disPosSub2 = BF_find_entry(cs,disPos2);
    tmp = zeros(Nx,numel(cs));
    tmp(cs+((1:numel(cs))-1)*Nx) = 1;
    data = fun(tmp);
    data = BF_ModeDecom(data,opt);
    M1alla = abs(data);
    data = data./abs(data);
    M1all = imag(log(data));
    
    Ridxall = []; Cidxall = []; Ridxalla = []; Cidxalla = [];
    for cntc = 1:opt.numCom
        M1 = M1all(:,:,cntc); M2 = M2all(:,:,cntc);
        [M1,M2] = BF_Phase_Correction_Mat_1D(M1,M2,cs,rs,disPos1,disPos2,disPosSub1,disPosSub2,2*pi);
        
        [~,R2,E2] = qr(M2',0);
        Cidx{cntc} = E2(find(abs(diag(R2))>tol*abs(R2(1)))<=tR);
        [~,R1,E1] = qr(M1',0);
        Ridx{cntc} = E1(find(abs(diag(R1))>tol*abs(R1(1)))<=tR);
        Ridxall = [Ridxall Ridx{cntc}]; Cidxall = [Cidxall Cidx{cntc}];
        
        M1 = M1alla(:,:,cntc); M2 = M2alla(:,:,cntc);
        [~,R2,E2] = qr(M2',0);
        Cidxa{cntc} = E2(find(abs(diag(R2))>tol*abs(R2(1)))<=tR);
        [~,R1,E1] = qr(M1',0);
        Ridxa{cntc} = E1(find(abs(diag(R1))>tol*abs(R1(1)))<=tR);
        Ridxalla = [Ridxalla Ridxa{cntc}]; Cidxalla = [Cidxalla Cidxa{cntc}];
    end
    
    if 0
        %rows again
        rs = BF_RandSample(Nx,tR);
        rs = unique([[disPos1,disPos1+1,disPos1+2,rs] Ridxall Ridxalla]);
        disPosSub1 = BF_find_entry(rs,disPos1);
        tmp = zeros(Np,numel(rs));
        tmp(rs+((1:numel(rs))-1)*Np) = 1;
        data = funt(tmp);
        data = BF_ModeDecom(data,opt);
        M2alla = abs(data);
        if type == 0
            data = data./abs(data);
        else
            data = conj(data)./abs(data);
        end
        M2all = imag(log(data));
        
        %columns
        cs = BF_RandSample(Np,tR);
        cs = unique([[disPos2,disPos2+1,disPos2+2,cs] Cidxall Cidxalla]);
        disPosSub2 = BF_find_entry(cs,disPos2);
        tmp = zeros(Nx,numel(cs));
        tmp(cs+((1:numel(cs))-1)*Nx) = 1;
        data = fun(tmp);
        data = BF_ModeDecom(data,opt);
        M1alla = abs(data);
        data = data./abs(data);
        M1all = imag(log(data));
        
        Ridxall = []; Cidxall = []; Ridxalla = []; Cidxalla = [];
        for cntc = 1:opt.numCom
            M1 = M1all(:,:,cntc); M2 = M2all(:,:,cntc);
            [M1,M2] = BF_Phase_Correction_Mat_1D(M1,M2,cs,rs,disPos1,disPos2,disPosSub1,disPosSub2,2*pi);
            
            [~,R2,E2] = qr(M2',0);
            Cidx{cntc} = E2(find(abs(diag(R2))>tol*abs(R2(1)))<=tR);
            [~,R1,E1] = qr(M1',0);
            Ridx{cntc} = E1(find(abs(diag(R1))>tol*abs(R1(1)))<=tR);
            Ridxall = [Ridxall Ridx{cntc}]; Cidxall = [Cidxall Cidx{cntc}];
            
            M1 = M1alla(:,:,cntc); M2 = M2alla(:,:,cntc);
            [~,R2,E2] = qr(M2',0);
            Cidxa{cntc} = E2(find(abs(diag(R2))>tol*abs(R2(1)))<=tR);
            [~,R1,E1] = qr(M1',0);
            Ridxa{cntc} = E1(find(abs(diag(R1))>tol*abs(R1(1)))<=tR);
            Ridxalla = [Ridxalla Ridxa{cntc}]; Cidxalla = [Cidxalla Cidxa{cntc}];
        end
    end
else
    for cntc = 1:opt.numCom
        Ridx{cntc} = 1:Nx;
        Cidx{cntc} = 1:Np;
        Ridxa{cntc} = 1:Nx;
        Cidxa{cntc} = 1:Np;
    end
    Ridxall = 1:Nx;
    Cidxall = 1:Np;
    Ridxalla = 1:Nx;
    Cidxalla = 1:Np;
end

%rows
Ridxall = [disPos1,disPos1+1,disPos1+2, Ridxall, Ridxalla]; Ridxall = unique(Ridxall);
disPosSub1 = BF_find_entry(Ridxall,disPos1);
tmp = zeros(Np,numel(Ridxall));
tmp(Ridxall+((1:numel(Ridxall))-1)*Np) = 1;
data = funt(tmp);
data = BF_ModeDecom(data,opt);
MRalla = abs(data);
if type == 0
    data = data./abs(data);
else
    data = conj(data)./abs(data);
end
MRall = imag(log(data));

%columns
Cidxall = [disPos2,disPos2+1,disPos2+2, Cidxall, Cidxalla]; Cidxall = unique(Cidxall);
disPosSub2 = BF_find_entry(Cidxall,disPos2);
tmp = zeros(Nx,numel(Cidxall));
tmp(Cidxall+((1:numel(Cidxall))-1)*Nx) = 1;
data = fun(tmp);
data = BF_ModeDecom(data,opt);
MCalla = abs(data);
data = data./abs(data);
MCall = imag(log(data));

M1 = cell(opt.numCom,1); M2 = cell(opt.numCom,1);
QC = cell(opt.numCom,1); QR = cell(opt.numCom,1);
M1a = cell(opt.numCom,1); M2a = cell(opt.numCom,1);
QCa = cell(opt.numCom,1); QRa = cell(opt.numCom,1);
for cntc = 1:opt.numCom
    MC = MCall(:,:,cntc); MR = MRall(:,:,cntc);
    [MC,MR] = BF_Phase_Correction_Mat_1D(MC,MR,Cidxall,Ridxall,disPos1,disPos2,disPosSub1,disPosSub2,2*pi);
    [QC{cntc},~,~] = qr(MC,0);
    [QR{cntc},~,~] = qr(MR,0);
    
    MC = MCalla(:,:,cntc); MR = MRalla(:,:,cntc);
    [QCa{cntc},~,~] = qr(MC,0);
    [QRa{cntc},~,~] = qr(MR,0);
end
if(tR<Np && tR<Nx)
    cs = BF_RandSample(Np,tR);
    cs = unique([cs Cidxall]); cs = sort(cs);
    rs = BF_RandSample(Nx,tR);
    rs = unique([rs Ridxall]); rs = sort(rs);
else
    cs = 1:Np;
    rs = 1:Nx;
end
disPosSub2 = BF_find_entry(cs,disPos2);
disPosSub1 = BF_find_entry(rs,disPos1);

for cntc = 1:opt.numCom
    M1{cntc} = QC{cntc}(rs,:);
    M2{cntc} = QR{cntc}(cs,:);
    M1a{cntc} = QCa{cntc}(rs,:);
    M2a{cntc} = QRa{cntc}(cs,:);
end

tmp = zeros(Nx,numel(cs));
tmp(cs+((1:numel(cs))-1)*Nx) = 1;
data = fun(tmp);
data = BF_ModeDecom(data,opt);
M1salla = abs(data);
data = data./abs(data);
M1sall = imag(log(data));

tmp = zeros(Np,numel(rs));
tmp(rs+((1:numel(rs))-1)*Np) = 1;
data = funt(tmp);
data = BF_ModeDecom(data,opt);
if type == 0
    data = data./abs(data);
else
    data = conj(data)./abs(data);
end
M2sall = imag(log(data));

U = cell(opt.numCom,1); V = cell(opt.numCom,1); S = cell(opt.numCom,1);
Ua = cell(opt.numCom,1); Va = cell(opt.numCom,1); Sa = cell(opt.numCom,1);
for cntc = 1:opt.numCom
    M1s = M1sall(:,:,cntc); M2s = M2sall(:,:,cntc);
    [M1s,~] = BF_Phase_Correction_Mat_1D(M1s,M2s,cs,rs,disPos1,disPos2,disPosSub1,disPosSub2,2*pi);
    M3 = M1s(rs,:);
    
    MD = pinv(M1{cntc}) * (M3* pinv(M2{cntc}'));
    [U{cntc},S{cntc},V{cntc}] = BF_svdtrunc_rank(MD,mR,tol);
    U{cntc} = QC{cntc}*U{cntc};
    V{cntc} = QR{cntc}*V{cntc};
    
    M1s = M1salla(:,:,cntc); 
    M3 = M1s(rs,:);
    
    MD = pinv(M1a{cntc}) * (M3* pinv(M2a{cntc}'));
    [Ua{cntc},Sa{cntc},Va{cntc}] = BF_svdtrunc_rank(MD,mR,tol);
    Ua{cntc} = QCa{cntc}*Ua{cntc};
    Va{cntc} = QRa{cntc}*Va{cntc};
end
end
