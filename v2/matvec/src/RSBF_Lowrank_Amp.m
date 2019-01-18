function [U,S,V] = RSBF_Lowrank_Amp(Nx,Np,fun,funt,tol,tR,mR,dim,isAL)
%
% Input:
% Suppose B is a complementary low-rank matrix
% fun is the function handle s.t. fun(idx) = B(:,idx)
% funt is the function handle s.t. funt(idx) = B(idx,:).'
% [Nx,Np] is the size of the matrix fun
% tol is the acccuracy tollerence
% tR is the test rank
% mR is the maximum rank, tR>=mR
% dim is the dimension of the problem
% isAL, if the input function handle is real, whether we use the analytic
% part of the data.
%
% Output:
% Low-rank factorization s.t. U*S*V' \approx abs(B)
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

if(tR<Np && tR<Nx)
    %rows
    rs = BF_RandSample(Nx,tR);
    data = funt(rs);
    if isAL
        data = BF_analytic(data,2);
    end
    data = abs(data);
    M2 = data;
    
    %columns
    cs = BF_RandSample(Np,tR);
    data = fun(cs);
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
    rs = unique([rs Ridx]);
    data = funt(rs);
    if isAL
        data = BF_analytic(data,2);
    end
    data = abs(data);
    M2 = data;
    
    %columns
    cs = BF_RandSample(Np,tR);
    cs = unique([cs Cidx]);
    data = fun(cs);
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
data = funt(Ridx);
if isAL
    data = BF_analytic(data,2);
end
data = abs(data);
MR = data;

%columns
data = fun(Cidx);
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
    cs = unique([cs Cidx]); cs = sort(cs);
    rs = BF_RandSample(Nx,tR);
    rs = unique([rs Ridx]); rs = sort(rs);
else
    cs = 1:Np;
    rs = 1:Nx;
end

M1 = QC(rs,:);
M2 = QR(cs,:);
data = fun(cs);
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
