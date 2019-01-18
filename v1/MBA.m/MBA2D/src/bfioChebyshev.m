function u = bfioChebyshev(Nx,Nk,SL,EL,EPS,fun,f,mats,dir,dirlev,stoplev,lev)
if Nk<=2^stoplev
    % Direct evaluation
    % When N<64, bfio is not acurate
    kg = -Nk/2:Nk/2-1;
    [k1,k2] = ndgrid(kg);
    ks = [k1(:)'; k2(:)'];
    xg = [0:1/Nx:(Nx-1)/Nx];
    [x1,x2] = ndgrid(xg);
    xs = [x1(:)'; x2(:)'];
    u = zeros(Nx^2,1);
    for cnt = 1:Nk^2
        u = u + fun(xs,ks(:,cnt))*(f(cnt));
    end
    u = reshape(u,Nx,Nx);
else
    u = bfio_eval(Nx,Nk,SL,EL,EPS,fun,f,mats,dir,dirlev{lev});
    tmp = bfioChebyshev(Nx,Nk/2,SL-1,EL,EPS,fun,f(end/4+1:3*end/4,end/4+1:3*end/4),mats,dir,dirlev,stoplev,lev+1);
    u = u + tmp;
end