function fun = directjac2d(nts,s,t,wghts,wghtt,n,da,db,c,flag)
if nargin < 10
    flag = [];
end
if nts < 2^12
   it = 10;
else
   it = 28;
end
vals = jacrecur(nts,s,wghts,it-1,da,db);
valt = jacrecur(nts,t,wghtt,it-1,da,db);
m = length(n);
nu = [it:nts-1]';
nt = zeros(nts,1);
if  flag ~= 'FULL'
    fun = @(c)dir2d1(c);

else
    Js = interpjac1(nt,s,nu,da,db,-1);
    Js = real(Js.*exp(1i*s*nu'));
    Js = [vals Js].';
    Jt = interpjac1(nt,t,nu,da,db,-1);
    Jt = real(Jt.*exp(1i*t*nu'));
    Jt = [valt Jt];
    fun = @(c)dir2d2(c);
end


    function y = dir2d1(c)
        c = c(:);
        y = zeros(m,1);
        for i = 1:m
            jj = floor((n(i)-0.5)/nts);
            ii = n(i) - nts*jj;
            x = interpjac1(nt,s(jj+1),nu,da,db,-1);
            %x1 = exp(1i*2*pi/nts*floor(s(jj+1)*nts/2/pi)*nu');
            x1 = exp(1i*s(jj+1)*nu');
            x = real(x.*x1);
            z = interpjac1(nt,t(ii),nu,da,db,-1);
            %z1 = exp(1i*2*pi/nts*floor(t(ii)*nts/2/pi)*nu');
            z1 = exp(1i*t(ii)*nu');
            z = real(z.*z1);
            x = [vals(jj+1,:) x];
            z = [valt(ii,:) z];
            y(i) = real(kron(x,z)*c);

        end
    end

    function y = dir2d2(c)
        y = Jt*c*Js;
        y = y(:);
    end
end
