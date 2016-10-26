function res = funF(N,x,k)
  phs = 2*pi*(x(1,:).'*k(1,:)+x(2,:).'*k(2,:));
  res = cos(phs) + i*sin(phs);
  
  
  
  