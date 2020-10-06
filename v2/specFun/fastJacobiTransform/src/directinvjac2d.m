function fun = directinvjac2d(nts,t,s,wghtt,wghts,n,da,db,c,flag)
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
    fun = @(c)dirinv2d1(c);

else
    WJs = interpjac1(nt,s,nu,da,db,-1);
    WJs = real(WJs.*exp(1i*s*nu'));
    WJs = diag(sqrt(wghts))*[vals WJs];
    WJs = WJs.';
    WJt = interpjac1(nt,t,nu,da,db,-1);
    WJt = real(WJt.*exp(1i*t*nu'));
    WJt = diag(sqrt(wghtt))*[valt WJt];
    fun = @(c)dirinv2d2(c);

end

    function y = dirinv2d1(c)
        y = zeros(m,1);
        c = c(:);
        c = kron(sqrt(wghtt),sqrt(wghts)).*c;

        for i = 1:m
            jj = floor((n(i)-0.5)/nts);
            ii = n(i) - nts*jj;
            if  jj < it
                x = valt(:,jj+1);
                x = sqrt(wghtt).*x;
                x = x.';
            else
                x = interpjac1(nt,t,jj,da,db,-1);
                x = real(sqrt(wghtt).*x.*exp(1i*t*jj));
                x = x.';
            end
            if  ii < it+1
                z = vals(:,ii);
                z = sqrt(wghts).*z;
                z = z.';
            else
                z = interpjac1(nt,s,ii-1,da,db,-1);
                z = real(sqrt(wghts).*z.*exp(1i*s*(ii-1)));
                z = z.';
            end
            y(i) = real(kron(x,z)*c);
        end
    end
    function y = dirinv2d2(c)
        c = diag(sqrt(wghts))*c*diag(sqrt(wghtt));
        y = WJs*c*WJt;
        y = y(:);
    end
end
