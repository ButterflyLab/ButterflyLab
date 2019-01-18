function y = polyFuncLege(x,k)
[valp,valq,alpha,alphader,vallogp,vallogq] = fastALegendre(k,0,x);
y = exp(1i*alpha);
end
