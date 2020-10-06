function y = directinvjac1d(nts,t,wghtt,n,da,db,c)
if nts < 2^12
   it = 10;
else
   it = 28;
end
valt = jacrecur(nts,t,wghtt,it-1,da,db);
m = length(n);
y = zeros(m,size(c,2));
nu = [it:nts-1]';
nt = zeros(nts,1);
c = sqrt(wghtt).*c;

for i = 1:m
    if  n(i) <= it
        v = (sqrt(wghtt).*valt(:,n(i))).';
        y(i,:) = v*c;
    else
        x = interpjac1(nt,t,n(i)-1,da,db,-1);
        x = real(sqrt(wghtt).*x.*exp(1i*t*(n(i)-1)));
        x = x.';
        y(i,:) = x*c;
    end
end
end