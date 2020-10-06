function fun = directinvjac3d(nts,r,s,t,wghtr,wghts,wghtt,n,da,db,c,flag)
if nargin < 10
    flag = [];
end

if nts < 2^12
   it = 10;
else
   it = 28;
end
valr = jacrecur(nts,r,wghtr,it-1,da,db);
vals = jacrecur(nts,s,wghts,it-1,da,db);
valt = jacrecur(nts,t,wghtt,it-1,da,db);
m = length(n);
nu = [it:nts-1]';
nt = zeros(nts,1);
if  flag == 'FULL'
    fun = @(c)dirinv3d1(c);

else
    WJr = interpjac1(nt,r,nu,da,db,-1);
    WJr = real(WJr.*exp(1i*r*nu'));
    WJr = diag(sqrt(wghtr))*[valr WJr];
    WJr = WJr.';
    WJs = interpjac1(nt,s,nu,da,db,-1);
    WJs = real(WJs.*exp(1i*s*nu'));
    WJs = diag(sqrt(wghts))*[vals WJs];
    WJs = WJs.';
    WJt = interpjac1(nt,t,nu,da,db,-1);
    WJt = real(WJt.*exp(1i*t*nu'));
    WJt = diag(sqrt(wghtt))*[valt WJt];
    WST = sparse(diag(kron(sqrt(wghts),sqrt(wghtt))));
    fun = @(c)dirinv3d2(c);

end
    function w = dirinv3d1(c)
        w = zeros(m,1);
        c = c(:);
        c = kron(kron(sqrt(wghtr),sqrt(wghts)),sqrt(wghtt)).*c;
        for i = 1:m
            kk = floor((n(i)-0.5)/nts/nts);
            jj = floor((n(i)-kk*nts*nts-0.5)/nts);
            ii = n(i)-kk*nts*nts-nts*jj;
            if  kk < it
                x = valr(:,kk+1);
                x = sqrt(wghtr).*x;
                x = x.';
            else
                x = interpjac1(nt,r,kk,da,db,-1);
                x = real(sqrt(wghtr).*x.*exp(1i*r*kk));
                x = x.';
            end
            if  jj < it
                y = vals(:,jj+1);
                y = sqrt(wghts).*y;
                y = y.';
            else
                y = interpjac1(nt,s,jj,da,db,-1);
                y = real(sqrt(wghts).*y.*exp(1i*s*jj));
                y = y.';
            end
            if  ii < it+1
                z = valt(:,ii);
                z = sqrt(wghtt).*z;
                z = z.';
            else
                z = interpjac1(nt,t,ii-1,da,db,-1);
                z = real(sqrt(wghtt).*z.*exp(1i*t*(ii-1)));
                z = z.';
            end
            w(i) = real(kron(kron(x,y),z)*c);
        end
    end

    function y = dirinv2d(c)
        y = WJt*c*WJs;
        y = y(:);
    end
    function y = dirinv3d2(c)
         c = reshape(c,nts^2,nts);
         c = WST*c*sparse(diag(sqrt(wghtr)));
         c = reshape(nts,nts,nts);
         y = zeros(nts^2,nts);
         for i = 1:nts
             y(:,i) = dirinv2d(c(:,:,i));
         end
         y = y*WJr;
         y = y(:);
    end
end
