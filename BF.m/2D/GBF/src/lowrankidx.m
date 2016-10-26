function V = lowrankidx(x,p,fun,Ridx,Vtmp)

Nx = size(x,1);
Np = size(p,1);

if(Nx==0 || Np==0)
    V = zeros(Np,0);
    return;
end

%get rows
MR = fun(x(Ridx,:),p);

%get middle matrix
[QR,~,~] = qr(MR',0);

V = QR*Vtmp;

end