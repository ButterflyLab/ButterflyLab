function res = aun0(N,x,k)
  nx = size(x,2);
  nk = size(k,2);
  %-----------------
  %xk = (x(1,:)'*k(1,:) + x(2,:)'*k(2,:));
  
  kr = sqrt(k(1,:).^2 + k(2,:).^2);
  sx = (2 + sin(2*pi*x(1,:)).*sin(2*pi*x(2,:)))/3; %CHECK C++ TOO
  
  bad = find(kr<eps); %zero entries
  arg = (2*pi)*sx(:)*kr;
  res = besselh(0,arg);
  res(:,bad) = 0;
  
  res = res .* exp(-i*arg);
end