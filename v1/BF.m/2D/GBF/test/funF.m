function res = funF(x,k)

xk = (x(:,1)*k(:,1)' + x(:,2)*k(:,2)');

tmp = (2*pi)* (xk);

res = exp(1i*tmp);

end
