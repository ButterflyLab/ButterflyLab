function res = fun0_2D(x,k)

xk = (x(:,1)*k(:,1)' + x(:,2)*k(:,2)');
sx = (2 + sin(2*pi*x(:,1)).*sin(2*pi*x(:,2)))/32;
cx = (2 + cos(2*pi*x(:,1)).*cos(2*pi*x(:,2)))/32;
kr = sqrt(sx.^2*(k(:,1).^2)' + cx.^2*(k(:,2).^2)');

tmp = (2*pi)* (xk + kr);

res = complex(cos(tmp),sin(tmp));

end
