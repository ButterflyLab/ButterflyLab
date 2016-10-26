function res = funF(N,x,k)
  phs = 2*pi*(x(1,:).'*k(1,:)+x(2,:).'*k(2,:)+x(3,:).'*k(3,:));
  res = cos(phs) + i*sin(phs);
  
  
  
  