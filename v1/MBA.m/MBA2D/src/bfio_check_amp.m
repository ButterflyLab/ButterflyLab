function relerr = bfio_check_amp(N,fun,aun,f,u,NC)
    kg = [-N/2:N/2-1];
    [k1,k2] = ndgrid(kg);
    ks = [k1(:)'; k2(:)'];
    app = zeros(NC,1);
    ext = zeros(NC,1);
    for g=1:NC
        x1 = floor(rand(1)*N);
        x2 = floor(rand(1)*N);
        app(g) = u(x1+1,x2+1);
        trg = [x1/N; x2/N];
        tmp = aun(trg,ks).*fun(trg,ks);
        ext(g) = sum(tmp*f(:));
    end
    err = app-ext;
    relerr = norm(err)/norm(ext);
end