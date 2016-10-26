function relerr = fbf_check(N,fun,f,xx,kk,u,NC)

app = zeros(NC,1);
ext = zeros(NC,1);
for g=1:NC
    x = floor(rand(1)*N+1);
    app(g) = u(x);
    ext(g) = sum(fun(xx(x,:),kk)*(f(:)));
end
err = app-ext;
relerr = norm(err)/norm(ext);

end