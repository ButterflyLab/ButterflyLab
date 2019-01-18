function u = bfioChebyshevAmp(Nx,Nk,SL,EL,EPS,fun,aun,arp,f,mats,dir,dirlev,stoplev,lev)
    psidx = arp{1};
    ptidx = arp{2};
    mid = arp{3};

    [k1s,k2s] = ndgrid([-Nk/2:Nk/2-1]);
    k1s = k1s(:)';  k2s = k2s(:)';
    [x1s,x2s] = ndgrid([0:Nx-1]/Nx);
    x1s = x1s(:)';  x2s = x2s(:)';
    pt = [x1s; x2s];
    ps = [k1s; k2s];

    M2 = mid * aun(ptidx,ps);
    M1 = aun(pt,psidx);

    u = zeros(Nx,Nx);
    for k=1:size(M1,2)
        tmpf = reshape(M2(k,:),Nk,Nk).*f;
        tmpu = bfioChebyshev(Nx,Nk,SL,EL,EPS,fun,tmpf,mats,dir,dirlev,stoplev,lev);
        u = u + reshape(M1(:,k),Nx,Nx).*tmpu;
    end
end