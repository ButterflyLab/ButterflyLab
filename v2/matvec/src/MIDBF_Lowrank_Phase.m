function [U,S,V,Dr,Dc] = MIDBF_Lowrank_Phase(Nx,Np,fun,funt,tol,tR,mR,thre,P1,P2)

% Input:
% Suppose Phi is a complementary low-rank matrix
% [Nx,Np] is the size of the matrix fun
% fun is the function handle s.t. fun(idx) = Phi(:,idx)
% funt is the function handle s.t. funt(idx) = Phi(idx,:).'
% tol is the acccuracy tollerence
% tR is the test rank
% mR is the maximum rank, tR>=mR
% thre is the threshold for discontinuity detection
% P1 and P2 are recovery path matrices by row and by column
%
% Output:
% Low-rank factorization s.t. exp(1i*U*S*V') \approx Phi
% 
% O(N log N) operation and memory complexity.
% Reference: Multidimensional Phase Recovery and Interpolative 
% Decomposition Butterfly Factorization, preprint, 2019.
% Copyright 2018 by Ze Chen, Juan Zhang, Ken Ho, and Haizhao Yang

iter = 2; % randomly sampling iteration

if(Nx==0 || Np==0)
    U = zeros(Nx,0);
    S = zeros(0,0);
    V = zeros(Np,0);
    return;
end

% test discontinuity along row
rs = BF_RandSample(Nx, tR);
dis2 = 1;
P2c = P2;
for cnt = 1:tR
    data = funt(rs(cnt));
    M2 = imag(log(data));
    der2nd = mod(M2(P2c(:,1)) - M2(P2c(:,2)) + pi, 2*pi) - pi;
    dis = P2c(abs(der2nd)>thre, 2);
    P2c(abs(der2nd)>thre, :) = [];
    dis = dis(:)';
    dis2 = unique([dis2,dis]);
end
Dc = length(dis2) - 1;

% test discontinuity along column
cs = BF_RandSample(Np, tR);
dis1 = 1;
P1c = P1;
for cnt = 1:tR
    data = fun(cs(cnt));
    M1 = imag(log(data));
    der2nd = mod(M1(P1c(:,1)) - M1(P1c(:,2)) + pi, 2*pi) - pi;
    dis = P1c(abs(der2nd)>thre, 2);
    P1c(abs(der2nd)>thre, :) = [];
    dis = dis(:)';
    dis1 = unique([dis1,dis]);
end
Dr = length(dis1) - 1;

if(tR<Np && tR<Nx)
    % rows
    rs = unique([1, dis1, BF_RandSample(Nx-1, tR-1)+1]);
    M2 = imag(log(funt(rs)));
    % columns
    cs = unique([1, dis2, BF_RandSample(Np-1, tR-1)+1]);
    M1 = imag(log(fun(cs)));
    
    [M1,M2] = BF_Phase_Correction_Mat_m(M1,M2,cs,rs,P1,P1c,P2,P2c,2*pi);
    
    [~,R2,E2] = qr(M2',0);
    Cidx = E2(find(abs(diag(R2))>tol*abs(R2(1)))<=tR);
    [~,R1,E1] = qr(M1',0);
    Ridx = E1(find(abs(diag(R1))>tol*abs(R1(1)))<=tR);
    
    for i = 1 : iter - 1
        % rows again
        rs = unique([1, dis1, Ridx, BF_RandSample(Nx-1, tR-1)+1]);
        M2 = imag(log(funt(rs)));
        % columns again
        cs = unique([1, dis2, Cidx, BF_RandSample(Np-1, tR-1)+1]);
        M1 = imag(log(fun(cs)));

        [M1,M2] = BF_Phase_Correction_Mat_m(M1,M2,cs,rs,P1,P1c,P2,P2c,2*pi);

        [~,R2,E2] = qr(M2',0);
        Cidx = E2(find(abs(diag(R2))>tol*abs(R2(1)))<=tR);
        [~,R1,E1] = qr(M1',0);
        Ridx = E1(find(abs(diag(R1))>tol*abs(R1(1)))<=tR);
    end
else
    Ridx = 1:Nx;
    Cidx = 1:Np;
end

% rows
Ridx = unique([1, dis1, Ridx]);
MR = imag(log(funt(Ridx)));
% columns
Cidx = unique([1, dis2, Cidx]);
MC = imag(log(fun(Cidx)));

[MC,MR] = BF_Phase_Correction_Mat_m(MC,MR,Cidx,Ridx,P1,P1c,P2,P2c,2*pi);

% get middle matrix
[QC,~,~] = qr(MC,0);
[QR,~,~] = qr(MR,0);

if(tR<Np && tR<Nx)
    rs = unique([1, dis1, Ridx, BF_RandSample(Nx-1, tR-1)+1]);
    cs = unique([1, dis2, Cidx, BF_RandSample(Np-1, tR-1)+1]);
else
    cs = 1:Np;
    rs = 1:Nx;
end

M1 = QC(rs,:);
M2 = QR(cs,:);
M1s = imag(log(fun(cs)));
M2s = imag(log(funt(rs)));

[M1s,~] = BF_Phase_Correction_Mat_m(M1s,M2s,cs,rs,P1,P1c,P2,P2c,2*pi);
M3 = M1s(rs,:);

MD = pinv(M1) * (M3 * pinv(M2'));
[U,S,V] = BF_svdtrunc_rank(MD,mR,tol);
U = QC*U;
V = QR*V;

r = find(abs(diag(S))>tol*abs(S(1)))<=mR;
S = S(r, r);
U = U(:, r);
V = V(:, r);

end
