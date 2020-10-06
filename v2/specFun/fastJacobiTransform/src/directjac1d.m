function y = directjac1d(nts,s,wghts,n,da,db,c)
if nts < 2^12
   it = 10;
else
   it = 28;
end
vals = jacrecur(nts,s,wghts,it-1,da,db);
m = length(n);
y = zeros(m,size(c,2));
nu = [it:nts-1]';
nt = zeros(nts,1);

for i = 1:m
    x = interpjac1(nt,s(n(i)),nu,da,db,-1);
    %x1 = exp(1i*2*pi/nts*floor(s(jj+1)*nts/2/pi)*nu');
    x1 = exp(1i*s(n(i))*nu');
    x = real(x.*x1);
    x = [vals(n(i),:) x];
    y(i,:) = x*c;
end
end