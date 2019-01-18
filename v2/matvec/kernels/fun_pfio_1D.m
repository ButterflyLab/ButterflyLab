function tmp = fun_pfio_1D(x,k)

xk = x*k';
sx = (2 + sin(2*pi*x))/8;
tmp = (2*pi)* (xk + sx*abs(k'));

end
