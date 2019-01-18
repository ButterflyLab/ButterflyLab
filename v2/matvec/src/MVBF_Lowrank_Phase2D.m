function [U,S,V] = MVBF_Lowrank_Phase(Nx,Np,fun,funt,tol,tR,mR,dim,thre,type,isAL)
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
% isAL, if the input function handle is real, whether we use the analytic
% part of the data.
%
% Output:
% Low-rank factorization s.t. exp(1i*U*S*V') \approx fun./abs(fun)
%
% O(N log N) operation and memory complexity.
% Reference: H. Yang, A Unified Framework for Oscillatory Integral 
% Transform: When to use NUFFT or Butterfly factorization? preprint, 2018.

if nargin < 11
    isAL = 0;
end
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
rs = BF_RandSample(Nx,tR);
disPos2 = [];
for cnt = 1:tR
    tmp = zeros(Nx,1); tmp(rs(cnt)) = 1;
    data = funt(tmp);
    if isAL
        data = BF_analytic(data,2);
    end
    if type == 0
        data = data./abs(data);
    else
        data = conj(data)./abs(data);
    end
    M2 = imag(log(data));
    M2 = BF_Phase_Correction_Vec_1(M2,1,2*pi,thre);
    der2nd = M2(1:end-2)+M2(3:end)-2*M2(2:end-1);
    dis = find(abs(der2nd)>thre)+2;
    dis = dis(:)';
    disPos2 = unique([disPos2,dis]);
end

% test discontinuity along column
cs = BF_RandSample(Np,tR);
disPos1 = [];
for cnt = 1:tR
    tmp = zeros(Np,1); tmp(cs(cnt)) = 1;
    data = fun(tmp);
    if isAL
        data = BF_analytic(data,2);
    end
    data = data./abs(data);
    M1 = imag(log(data));
    M1 = BF_Phase_Correction_Vec_1(M1,1,2*pi,thre);
    der2nd = M1(1:end-2)+M1(3:end)-2*M1(2:end-1);
    dis = find(abs(der2nd)>thre)+2;
    dis = dis(:)';
    disPos1 = unique([disPos1,dis]);
end

disPos1 = unique([1,disPos1]); 
disPos1 = setdiff(disPos1,disPos1+1);
disPos1 = setdiff(disPos1,disPos1+2);
disPos2 = unique([1,disPos2]); 
disPos2 = setdiff(disPos2,disPos2+1);
disPos2 = setdiff(disPos2,disPos2+2);

if(tR<Np && tR<Nx)
    %rows
    rs = BF_RandSample(Nx,tR);
    rs = [disPos1,disPos1+1,disPos1+2,rs]; rs = unique(rs); rs = sort(rs);
    disPosSub1 = BF_find_entry(rs,disPos1);
    tmp = zeros(Nx,numel(rs));
    tmp(rs+((1:numel(rs))-1)*Nx) = 1;
    data = funt(tmp);
    if isAL
        data = BF_analytic(data,2);
    end
    if type == 0
        data = data./abs(data);
    else
        data = conj(data)./abs(data);
    end
    M2 = imag(log(data));
    
    %columns
    cs = BF_RandSample(Np,tR);
    cs = [disPos2,disPos2+1,disPos2+2,cs]; cs = unique(cs); cs = sort(cs);
    disPosSub2 = BF_find_entry(cs,disPos2);
    tmp = zeros(Np,numel(cs));
    tmp(cs+((1:numel(cs))-1)*Np) = 1;
    data = fun(tmp);
    if isAL
        data = BF_analytic(data,2);
    end
    data = data./abs(data);
    M1 = imag(log(data));
    
    [M1,M2] = BF_Phase_Correction_Mat_1D(M1,M2,cs,rs,disPos1,disPos2,disPosSub1,disPosSub2,2*pi);
    
    [~,R2,E2] = qr(M2',0);
    Cidx = E2(find(abs(diag(R2))>tol*abs(R2(1)))<=tR);
    [~,R1,E1] = qr(M1',0);
    Ridx = E1(find(abs(diag(R1))>tol*abs(R1(1)))<=tR);
    
    %rows again
    rs = BF_RandSample(Nx,tR);
    rs = unique([[disPos1,disPos1+1,disPos1+2,rs] Ridx]);
    disPosSub1 = BF_find_entry(rs,disPos1);
    tmp = zeros(Nx,numel(rs));
    tmp(rs+((1:numel(rs))-1)*Nx) = 1;
    data = funt(tmp);
    if isAL
        data = BF_analytic(data,2);
    end
    if type == 0
        data = data./abs(data);
    else
        data = conj(data)./abs(data);
    end
    M2 = imag(log(data));
    
    %columns
    cs = BF_RandSample(Np,tR);
    cs = unique([[disPos2,disPos2+1,disPos2+2,cs] Cidx]);
    disPosSub2 = BF_find_entry(cs,disPos2);
    tmp = zeros(Np,numel(cs));
    tmp(cs+((1:numel(cs))-1)*Np) = 1;
    data = fun(tmp);
    if isAL
        data = BF_analytic(data,2);
    end
    data = data./abs(data);
    M1 = imag(log(data));
    
    [M1,M2] = BF_Phase_Correction_Mat_1D(M1,M2,cs,rs,disPos1,disPos2,disPosSub1,disPosSub2,2*pi);
    
    [~,R2,E2] = qr(M2',0);
    Cidx = E2(find(abs(diag(R2))>tol*abs(R2(1)))<=tR);
    [~,R1,E1] = qr(M1',0);
    Ridx = E1(find(abs(diag(R1))>tol*abs(R1(1)))<=tR);
else
    Ridx = 1:Nx;
    Cidx = 1:Np;
end

%rows
Ridx = [disPos1,disPos1+1,disPos1+2, Ridx]; Ridx = unique(Ridx);
disPosSub1 = BF_find_entry(Ridx,disPos1);
tmp = zeros(Nx,numel(Ridx));
tmp(Ridx+((1:numel(Ridx))-1)*Nx) = 1;
data = funt(tmp);
if isAL
    data = BF_analytic(data,2);
end
if type == 0
    data = data./abs(data);
else
    data = conj(data)./abs(data);
end
MR = imag(log(data));

%columns
Cidx = [disPos2,disPos2+1,disPos2+2, Cidx]; Cidx = unique(Cidx);
disPosSub2 = BF_find_entry(Cidx,disPos2);
tmp = zeros(Np,numel(Cidx));
tmp(Cidx+((1:numel(Cidx))-1)*Np) = 1;
data = fun(tmp);
if isAL
    data = BF_analytic(data,2);
end
data = data./abs(data);
MC = imag(log(data));

[MC,MR] = BF_Phase_Correction_Mat_1D(MC,MR,Cidx,Ridx,disPos1,disPos2,disPosSub1,disPosSub2,2*pi);

%get middle matrix
[QC,~,~] = qr(MC,0);
[QR,~,~] = qr(MR,0);

if(tR<Np && tR<Nx)
    cs = BF_RandSample(Np,tR);
    cs = unique([cs Cidx]); cs = sort(cs);
    rs = BF_RandSample(Nx,tR);
    rs = unique([rs Ridx]); rs = sort(rs);
else
    cs = 1:Np;
    rs = 1:Nx;
end
disPosSub2 = BF_find_entry(cs,disPos2);
disPosSub1 = BF_find_entry(rs,disPos1);

M1 = QC(rs,:);
M2 = QR(cs,:);
tmp = zeros(Np,numel(cs));
tmp(cs+((1:numel(cs))-1)*Np) = 1;
data = fun(tmp);
if isAL
    data = BF_analytic(data,2);
end
data = data./abs(data);
M1s = imag(log(data));

tmp = zeros(Nx,numel(rs));
tmp(rs+((1:numel(rs))-1)*Nx) = 1;
data = funt(tmp);
if isAL
    data = BF_analytic(data,2);
end
if type == 0
    data = data./abs(data);
else
    data = conj(data)./abs(data);
end
M2s = imag(log(data));

[M1s,~] = BF_Phase_Correction_Mat_1D(M1s,M2s,cs,rs,disPos1,disPos2,disPosSub1,disPosSub2,2*pi);
M3 = M1s(rs,:);

MD = pinv(M1) * (M3* pinv(M2'));
[U,S,V] = BF_svdtrunc_rank(MD,mR,tol);
U = QC*U;
V = QR*V;

end
