function res = fun_kernal_hel(k,x,y)

tmp = (2*pi) * (k*pdist2(x,y));
res = complex(cos(tmp),sin(tmp));

end