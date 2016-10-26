function res = aun1(x,k)

x1 = x(1,:);
x2 = x(2,:);
k1 = k(1,:);
k2 = k(2,:);

tmp1 = (2 + sin(2*pi*x1).*sin(2*pi*x2))/3;
tmp1 = tmp1';
tmp2 = (2 + cos(2*pi*x1).*cos(2*pi*x2))/3;
tmp2 = tmp2';

phs = sqrt( tmp1*k1.^2 + tmp2*k2.^2 );

bad = phs<eps;
arg = (2*pi)*phs;
res = besselh(0,arg);
res(bad) = 0;

res = res .* exp(-1i*arg);

end
