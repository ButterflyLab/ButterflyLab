function res = fun0(x,k)

xk = (x(:,1)*k(:,1)' + x(:,2)*k(:,2)');
sx = (2 + sin(2*pi*x(:,1)).*sin(2*pi*x(:,2)))/16;
cx = (2 + cos(2*pi*x(:,1)).*cos(2*pi*x(:,2)))/16;
kr = sqrt(sx.^2*(k(:,1).^2)' + cx.^2*(k(:,2).^2)');

tmp = (2*pi)* (xk + kr);

res = exp(1i*tmp);

end
