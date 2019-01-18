%function [U,S,V] = RSBF_Lowrank_Phase2D(Nx,Np,fun,funt,tol,tR,mR,thre,type,isAL,disPos1,disPos2)
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
% thre is the threshold for discontinuity detection
% type, when funt is transpose, type = 0; otherwise, type = 1.
% isAL, if the input function handle is real, whether we use the analytic
% part of the data.
%
% Output:
% Low-rank factorization s.t. exp(1i*U*S*V') \approx B./abs(B)
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

% % test discontinuity along row
% rs = BF_RandSample(Nx,tR);
% disPos2 = [];
% for cnt = 1:tR
%     data = funt(rs(cnt));
%     if isAL
%         data = BF_analytic(data,2);
%     end
%     if type == 0
%         data = data./abs(data);
%     else
%         data = conj(data)./abs(data);
%     end
%     M2 = imag(log(data));
%     M2 = BF_Phase_Correction_Vec_1(M2,1,2*pi,thre);
%     der2nd = M2(1:end-2)+M2(3:end)-2*M2(2:end-1);
%     dis = find(abs(der2nd)>thre)+2;
%     dis = dis(:)';
%     disPos2 = unique([disPos2,dis]);
% end
% 
% % test discontinuity along column
% cs = BF_RandSample(Np,tR);
% disPos1 = [];
% for cnt = 1:tR
%     data = fun(cs(cnt));
%     if isAL
%         data = BF_analytic(data,2);
%     end
%     data = data./abs(data);
%     M1 = imag(log(data));
%     M1 = BF_Phase_Correction_Vec_1(M1,1,2*pi,thre);
%     der2nd = M1(1:end-2)+M1(3:end)-2*M1(2:end-1);
%     dis = find(abs(der2nd)>thre)+2;
%     dis = dis(:)';
%     disPos1 = unique([disPos1,dis]);
% end
% 
% disPos1 = unique([1,disPos1]); 
% disPos1 = setdiff(disPos1,disPos1+1);
% disPos1 = setdiff(disPos1,disPos1+2);
% disPos2 = unique([1,disPos2]); 
% disPos2 = setdiff(disPos2,disPos2+1);
% disPos2 = setdiff(disPos2,disPos2+2);


% assume that we don't have discontinuity at all at this moment
if(tR<Np && tR<Nx)
    %rows
    rs = BF_RandSample(Nx,tR);
    rsi = reshape((repmat(disPos1(1:3,1),1,3)+repmat([0,1,2],3,1))',[9,1])';
    rs = [rsi,rs]; rs = unique(rs);
    outDisPos1 = find(disPos1(:,2)==1);
    disPosSub1 = BF_find_entry(rs,disPos1(outDisPos1,1));
    data = funt(rs);
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
    csi = reshape((repmat(disPos2(1:3,1),1,3)+repmat([0,1,2],3,1))',[9,1])';
    cs = [csi,cs]; cs = unique(cs);
    outDisPos2 = find(disPos2(:,2)==1);
    disPosSub2 = BF_find_entry(cs,disPos2(outDisPos2,1));
    data = fun(cs);
    if isAL
        data = BF_analytic(data,2);
    end
    data = data./abs(data);
    M1 = imag(log(data));
    
    [M1,M2] = BF_Phase_Correction_Mat_2D(M1,M2,cs,rs,disPos1,disPos2,disPosSub1,disPosSub2,2*pi);
    
    [~,R2,E2] = qr(M2',0);
    Cidx = E2(find(abs(diag(R2))>tol*abs(R2(1)))<=tR);
    [~,R1,E1] = qr(M1',0);
    Ridx = E1(find(abs(diag(R1))>tol*abs(R1(1)))<=tR);
    
    %rows again
    rs = BF_RandSample(Nx,tR);
    rs = unique([[disPos1,disPos1+1,disPos1+2,rs] Ridx]);
    disPosSub1 = BF_find_entry(rs,disPos1);
    data = funt(rs);
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
    data = fun(cs);
    if isAL
        data = BF_analytic(data,2);
    end
    data = data./abs(data);
    M1 = imag(log(data));
    
    [M1,M2] = BF_Phase_Correction_Mat_2D(M1,M2,cs,rs,disPos1,disPos2,disPosSub1,disPosSub2,2*pi);
    
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
data = funt(Ridx);
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
data = fun(Cidx);
if isAL
    data = BF_analytic(data,2);
end
data = data./abs(data);
MC = imag(log(data));

[MC,MR] = BF_Phase_Correction_Mat_2D(MC,MR,Cidx,Ridx,disPos1,disPos2,disPosSub1,disPosSub2,2*pi);

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
data = fun(cs);
if isAL
    data = BF_analytic(data,2);
end
data = data./abs(data);
M1s = imag(log(data));

data = funt(rs);
if isAL
    data = BF_analytic(data,2);
end
if type == 0
    data = data./abs(data);
else
    data = conj(data)./abs(data);
end
M2s = imag(log(data));

[M1s,~] = BF_Phase_Correction_Mat_2D(M1s,M2s,cs,rs,disPos1,disPos2,disPosSub1,disPosSub2,2*pi);
M3 = M1s(rs,:);

MD = pinv(M1) * (M3* pinv(M2'));
[U,S,V] = BF_svdtrunc_rank(MD,mR,tol);
U = QC*U;
V = QR*V;

end
