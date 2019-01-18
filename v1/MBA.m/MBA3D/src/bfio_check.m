function relerr = bfio_check(N,fun,f,u,NC)
    kg = [-N/2:N/2-1];
    [k1,k2,k3] = ndgrid(kg);
    ks = [k1(:)'; k2(:)'; k3(:)'];
    app = zeros(NC,1);
    ext = zeros(NC,1);
    for g=1:NC
        x1 = floor(rand(1)*N);
        x2 = floor(rand(1)*N);
        x3 = floor(rand(1)*N);
        app(g) = u(x1+1,x2+1,x3+1);
        ext(g) = sum(fun(N,[x1/N;x2/N;x3/N],ks)*(f(:)));
    end
    err = app-ext;
    relerr = norm(err)/norm(ext);
end