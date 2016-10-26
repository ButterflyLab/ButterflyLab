function [U,S,V,Ridx,Cidx,rs,cs] = lowrank(x,p,fun,tol,tR,mR)

Nx = size(x,1);
Np = size(p,1);

if(Nx==0 || Np==0)
    U = zeros(Nx,0);
    S = zeros(0,0);
    V = zeros(Np,0);
    Ridx = [];
    Cidx = [];
    rs = [];
    cs = [];
    return;
end

if(tR<Np && tR<Nx)
    %get columns
    rs = randsample(Nx,tR);
    M2 = fun(x(rs,:),p);
    [~,R2,E2] = qr(M2,0);
    Cidx = E2(find(abs(diag(R2))>tol*abs(R2(1)))<=tR);

    %get rows
    cs = randsample(Np,tR);
    cs = unique([cs' Cidx]);
    M1 = fun(x,p(cs,:));
    [~,R1,E1] = qr(M1',0);
    Ridx = E1(find(abs(diag(R1))>tol*abs(R1(1)))<=tR);

    %get columns again
    rs = randsample(Nx,tR);
    rs = unique([rs' Ridx]);
    M2 = fun(x(rs,:),p);
    [~,R2,E2] = qr(M2,0);
    Cidx = E2(find(abs(diag(R2))>tol*abs(R2(1)))<=tR);

    %get rows again
    cs = randsample(Np,tR);
    cs = unique([cs' Cidx]);
    M1 = fun(x,p(cs,:));
    [~,R1,E1] = qr(M1',0);
    Ridx = E1(find(abs(diag(R1))>tol*abs(R1(1)))<=tR);
else
    Ridx = 1:Nx;
    Cidx = 1:Np;
end

%get rows
MR = fun(x(Ridx,:),p);

%get columns
MC = fun(x,p(Cidx,:));

%get middle matrix
[QC,~,~] = qr(MC,0);
[QR,~,~] = qr(MR',0);

if(tR<Np && tR<Nx)
    cs = randsample(Np,tR);
    cs = unique([cs' Cidx]);
    rs = randsample(Nx,tR);
    rs = unique([rs' Ridx]);
else
    cs = 1:Np;
    rs = 1:Nx;
end

M1 = QC(rs,:);
M2 = QR(cs,:);
M3 = fun(x(rs,:),p(cs,:));
MD = pinv(M1) * (M3* pinv(M2'));
[U,S,V] = svdtrunc(MD,mR,tol);
U = QC*U;
V = QR*V;

end
