function relerr = pbf_explicit_check(N,fun,f,u,NC)

k = -N/2:N/2-1;
[k1,k2] = ndgrid(k);
k1 = k1(:);
k2 = k2(:);
kk = [k1 k2];
app = zeros(NC,1);
ext = zeros(NC,1);
for g=1:NC
    x1idx = floor(rand(1)*N+1);
    x2idx = floor(rand(1)*N+1);
    xxidx = x1idx + (x2idx-1)*N;
    x1 = (x1idx-1)/N;
    x2 = (x2idx-1)/N;
    xx = [x1 x2];
    app(g) = u(xxidx);
    ext(g) = sum(fun(xx,kk)*(f(:)));
end
err = app-ext;
relerr = norm(err)/norm(ext);

end
