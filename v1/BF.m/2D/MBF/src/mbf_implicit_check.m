function relerr = mbf_implicit_check(N,y,yy,NC)

app = zeros(NC,1);
ext = zeros(NC,1);
for g=1:NC
    x = floor(rand(1)*N*2+1);
    app(g) = yy(x);
    ext(g) = y(x);
end
err = app-ext;
relerr = norm(err)/norm(ext);

end