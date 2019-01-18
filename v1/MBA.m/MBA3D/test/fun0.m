function res = fun0(N,x,k)
  nx = size(x,2);
  nk = size(k,2);
  xk = (x(1,:)'*k(1,:) + x(2,:)'*k(2,:) + x(3,:)'*k(3,:));
  kr = sqrt(k(1,:).^2 + k(2,:).^2 + k(3,:).^2);
  sx = (2 + sin(2*pi*x(1,:)).*sin(2*pi*x(2,:)).*sin(2*pi*x(3,:)))/3;
  tmp = (2*pi)* (xk + sx(:)*kr);
  res = complex(cos(tmp),sin(tmp));
end