function res = fun0_1D(x,k)

xk = x*k';
sx = (2 + sin(2*pi*x))/8;
tmp = (2*pi)* (xk + sx*abs(k'));
res = complex(cos(tmp),sin(tmp));

end
