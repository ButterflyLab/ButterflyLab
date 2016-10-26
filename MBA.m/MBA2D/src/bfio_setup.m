function [arp,err] = bfio_setup(N,aun)
  
  [k1s,k2s] = ndgrid([-N/2:N/2-1]);
  k1s = k1s(:)';  k2s = k2s(:)';
  
  [x1s,x2s] = ndgrid([0:N-1]/N);
  x1s = x1s(:)';  x2s = x2s(:)';
  
  % sample
  pt = [x1s; x2s];
  ps = [k1s; k2s];
  
  m = size(pt,2);
  n = size(ps,2);
  
  EPS = 1e-6;
  NPK = 10;
  
  %step 1
  nr = NPK;
  rs = sort(ceil(rand(1,nr)*(m-nr))) + [1:nr];
  M1 = aun(pt(:,rs),ps);
  [Q1,R1,E1] = qr(M1,0); clear M1;
  gd = find(abs(diag(R1))>EPS*abs(R1(1))); 
  Q1 = Q1(:,gd);
  R1 = R1(gd,gd);
  idx1 = E1(gd);    psidx = ps(:,idx1); clear Q1 R1;
  
  %step 2
  nc = NPK;
  cs = sort(ceil(rand(1,nc)*(n-nc))) + [1:nc];  cs = unique([cs idx1]);
  M2 = aun(pt,ps(:,cs));
  M2 = M2';
  [Q2,R2,E2] = qr(M2,0); clear M2;
  gd = find(abs(diag(R2))>EPS*abs(R2(1)));
  Q2 = Q2(:,gd);
  R2 = R2(gd,gd);
  idx2 = E2(gd);  ptidx = pt(:,idx2); clear Q2 R2;
  
  nc = NPK;
  cs = sort(ceil(rand(1,nc)*(n-nc))) + [1:nc];
  cs = unique([idx1 cs]);
  nr = NPK;
  rs = sort(ceil(rand(1,nr)*(m-nr))) + [1:nr];
  rs = unique([idx2 rs]);
  
  M1 = aun(pt(:,rs),psidx);
  M2 = aun(ptidx,ps(:,cs));
  M3 = aun(pt(:,rs),ps(:,cs));
  
  eqn = psidx;
  chk = ptidx;
  mid = pinv(M1) * (M3 * pinv(M2));
  
  arp = {eqn chk mid};
  
  if(1) %check
    nc = NPK;
    cs = sort(ceil(rand(1,nc)*(n-nc))) + [1:nc];
    cs = unique([idx1 cs]);
    nr = NPK;
    rs = sort(ceil(rand(1,nr)*(m-nr))) + [1:nr];
    rs = unique([idx2 rs]);
    
    M1 = aun(pt(:,rs),psidx);
    M2 = aun(ptidx,ps(:,cs));
    M3 = aun(pt(:,rs),ps(:,cs));
    E = M3 - M1*mid*M2;
    err=norm(E)/norm(M3);
  end
