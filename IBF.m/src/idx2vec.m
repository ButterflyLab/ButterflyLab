function x = idx2vec(siz,iter)

n = length(siz);
x = zeros(1,n);

%assert(iter <= siz(end),'The iter is out of the range');

for i = 1:n-1
    x(i) = mod(iter-1,siz(i))+1;
    iter = floor((iter-0.5)/siz(i))+1;
end
x(n) = mod(iter-1,siz(n))+1;

end