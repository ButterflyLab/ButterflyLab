function [U,S,V] = BF_lowrankidx(x,p,fun,tol,mR,Ridx,Cidx,rs,cs)

Nx = size(x,1);
Np = size(p,1);

if(Nx==0 || Np==0)
    U = zeros(Nx,0);
    S = zeros(0,0);
    V = zeros(Np,0);
    return;
end

%get rows
MR = fun(x(Ridx,:),p);

%get columns
MC = fun(x,p(Cidx,:));

%get middle matrix
[QC,~,~] = qr(MC,0);
[QR,~,~] = qr(MR',0);

M1 = QC(rs,:);
M2 = QR(cs,:);
M3 = fun(x(rs,:),p(cs,:));
MD = pinv(M1) * (M3* pinv(M2'));
[U,S,V] = BF_svdtrunc_rank(MD,mR,tol);
U = QC*U;
V = QR*V;

end

