function res = funF(N,x,k)

x = (x-1)/N;
k = k-1-N/2;
tmp = (2*pi)* (x*k');
res = complex(cos(tmp),sin(tmp));

end
