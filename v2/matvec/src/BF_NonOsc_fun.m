function A = BF_NonOsc_fun(fun,x,k,x0,k0)
A = fun(x,k);
kk = repmat(k0,[size(k,1),1]);
A = A./fun(x,kk);
xx = repmat(x0,[size(x,1),1]);
A = A./fun(xx,k);
end