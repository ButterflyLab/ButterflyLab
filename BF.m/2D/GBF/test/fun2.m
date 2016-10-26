function res = fun2(x,k)

xk = (x(:,1)*k(:,1)' + x(:,2)*k(:,2)');
kr = sqrt(k(:,1).^2 + k(:,2).^2);
xr = sqrt(x(:,1).^2 + x(:,2).^2);

tmp = (2*pi)* (xk + xr*kr');

res = exp(1i*tmp);

end
