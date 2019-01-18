function idx = vec2idx(siz,vec)

%assert(all(vec <= siz),'The iter is out of the range');

n = length(siz);
siz = cumprod(siz);

idx = vec(1);
for i = 2:n
    idx = idx + (vec(i)-1)*siz(i-1);
end

end