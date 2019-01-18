function [U,S,V] = MVBF_Lowrank_Amp(Nx,Np,fun,funt,tol,tR,mR,dim,isAL)
%
% Input:
% fun is the compressed butterfly matrix in a form of a function handle
% funt is the transpose of fun in a form of a function handle
% [Nx,Np] is the size of the matrix fun
% tol is the acccuracy tollerence
% tR is the test rank
% mR is the maximum rank, tR>=mR
% dim is the dimension of the problem
% isAL, if the input function handle is real, whether we use the analytic
% part of the data.
%
% Output:
% Low-rank factorization s.t. U*S*V' \approx abs(fun)
%
% O(N log N) operation and memory complexity.
% Reference: H. Yang, A Unified Framework for Oscillatory Integral 
% Transform: When to use NUFFT or Butterfly factorization? preprint, 2018.

if nargin < 9, isAL = 0; end;
if nargin < 8, dim = 1; end;
if(Nx==0 || Np==0)
    U = zeros(Nx,0);
    S = zeros(0,0);
    V = zeros(Np,0);
    return;
end
mR = min(min(Nx,Np),mR);
tR = min(min(Nx,Np),tR);
grid = BF_Chey_grid(tR);
gd1 = round(grid*(Nx-min(Nx,tR)) + (0:min(Nx,tR)-1)')+1;
gd2 = round(grid*(Np-min(Np,tR)) + (0:min(Np,tR)-1)')+1;
disPos1m = gd1';
disPos2m = gd2';

if(tR<Np && tR<Nx)
    %rows
    rs = BF_RandSample(Nx,tR);
    rs = [rs,disPos1m]; rs = unique(rs); rs = sort(rs);
    tmp = zeros(Nx,numel(rs));
    tmp(rs+((1:numel(rs))-1)*Nx) = 1;
    data = funt(tmp);
    if isAL
        data = BF_analytic(data,2);
    end
    data = abs(data);
    M2 = data;
    
    %columns
    cs = BF_RandSample(Np,tR);
    cs = [cs,disPos2m]; cs = unique(cs); cs = sort(cs);
    tmp = zeros(Np,numel(cs));
    tmp(cs+((1:numel(cs))-1)*Np) = 1;
    data = fun(tmp);
    if isAL
        data = BF_analytic(data,2);
    end
    data = abs(data);
    M1 = data;
    
    [~,R2,E2] = qr(M2',0);
    Cidx = E2(find(abs(diag(R2))>tol*abs(R2(1)))<=tR);
    [~,R1,E1] = qr(M1',0);
    Ridx = E1(find(abs(diag(R1))>tol*abs(R1(1)))<=tR);
    
    %rows again
    rs = BF_RandSample(Nx,tR);
    rs = [rs,disPos1m]; rs = unique(rs); rs = sort(rs);
    rs = unique([rs Ridx]);
    tmp = zeros(Nx,numel(rs));
    tmp(rs+((1:numel(rs))-1)*Nx) = 1;
    data = funt(tmp);
    if isAL
        data = BF_analytic(data,2);
    end
    data = abs(data);
    M2 = data;
    
    %columns
    cs = BF_RandSample(Np,tR);
    cs = [cs,disPos2m]; cs = unique(cs); cs = sort(cs);
    cs = unique([cs Cidx]);
    tmp = zeros(Np,numel(cs));
    tmp(cs+((1:numel(cs))-1)*Np) = 1;
    data = fun(tmp);
    if isAL
        data = BF_analytic(data,2);
    end
    data = abs(data);
    M1 = data;
    
    [~,R2,E2] = qr(M2',0);
    Cidx = E2(find(abs(diag(R2))>tol*abs(R2(1)))<=tR);
    [~,R1,E1] = qr(M1',0);
    Ridx = E1(find(abs(diag(R1))>tol*abs(R1(1)))<=tR);
else
    Ridx = 1:Nx;
    Cidx = 1:Np;
end

%rows
Ridx = [Ridx,disPos1m]; Ridx = unique(Ridx);
tmp = zeros(Nx,numel(Ridx));
tmp(Ridx+((1:numel(Ridx))-1)*Nx) = 1;
data = funt(tmp);
if isAL
    data = BF_analytic(data,2);
end
data = abs(data);
MR = data;

%columns
Cidx = [Cidx,disPos2m]; Cidx = unique(Cidx);
tmp = zeros(Np,numel(Cidx));
tmp(Cidx+((1:numel(Cidx))-1)*Np) = 1;
data = fun(tmp);
if isAL
    data = BF_analytic(data,2);
end
data = abs(data);
MC = data;

%get middle matrix
[QC,~,~] = qr(MC,0);
[QR,~,~] = qr(MR,0);

if(tR<Np && tR<Nx)
    cs = BF_RandSample(Np,tR);
    cs = unique([cs,Cidx,disPos2m]); cs = sort(cs);
    rs = BF_RandSample(Nx,tR);
    rs = unique([rs,Ridx,disPos1m]); rs = sort(rs);
else
    cs = 1:Np;
    rs = 1:Nx;
end

M1 = QC(rs,:);
M2 = QR(cs,:);
tmp = zeros(Np,numel(cs));
tmp(cs+((1:numel(cs))-1)*Np) = 1;
data = fun(tmp);
if isAL
    data = BF_analytic(data,2);
end
data = abs(data);
M1s = data;

M3 = M1s(rs,:);

MD = pinv(M1) * (M3* pinv(M2'));
[U,S,V] = BF_svdtrunc_rank(MD,mR,tol);
U = QC*U;
V = QR*V;

end
