function res = bfio_prep_aux(gs,ts)
%gs - grid points
%ts - test points
  NG = size(gs,2);
  NT = size(ts,2);
  
  tmp = zeros(NG,NT);
  for a=1:NG
      gud = [1:a-1 a+1:NG];
      for b=1:NT
          cur = (ts(b)-gs(gud))./(gs(a)-gs(gud));
          tmp(a,b) = prod(cur);
      end
  end
  res = tmp;
  