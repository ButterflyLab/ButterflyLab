function res = fun1(N,x,k)
  nx = size(x,2);
  nk = size(k,2);
  x1 = x(1,:);  x2 = x(2,:);
  k1 = k(1,:);  k2 = k(2,:);
  xk = (x(1,:)'*k(1,:) + x(2,:)'*k(2,:));
  
  tmp1 = (2 + sin(2*pi*x1).*sin(2*pi*x2))/3; tmp1 = tmp1';
  tmp2 = (2 + cos(2*pi*x1).*cos(2*pi*x2))/3; tmp2 = tmp2';
  phs = sqrt( tmp1*k1.^2 + tmp2*k2.^2 );
  %res = exp( (2*pi*i) * (xk + phs));
  tmp = 2*pi*(xk+phs);  res = cos(tmp) + i*sin(tmp);
  
  %sx = (2 + sin(2*pi*x(1,:)/N).*sin(2*pi*x(2,:)/N));
  %res = exp( (2*pi*i)* (sx(:)*kr) );
  
  